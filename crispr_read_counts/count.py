import sys
from typing import Dict, Any
from . import error_msg, open_plain_or_gzipped_file, PLASMID_COUNT_HEADER
import re
import pysam


DNA_PATTERN = re.compile(r'^[ATGC]+$', flags=re.I)
SAFE_SEQ_FORMAT = re.compile(r'^F\d+$')
DUAL_LIBRARY_EXPECTED_HEADER = ['sgrna_left_id', 'sgrna_left_seq', 'sgrna_right_id', 'sgrna_right_seq', 'unique_id', 'gene_pair_id', 'target_id']

def count_single(args: Dict[str, Any]):

  # TODO: validate inputs before doing anything

  lib, targeted_genes = get_singe_guide_library(args['library'])

  plasmid, plas_name = None, None
  if args['plasmid']:
    plasmid, plas_name = get_plasmid_read_counts(args['plasmid'])

  sample, sample_name = get_sgrna_library_counts(args['dir'], args['ref'], lib, args['trim'], args['reverse_complement'])

  with open(args['output'], 'w') as f:
    if plas_name:
      f.write('\t'.join(['sgRNA', 'gene', f'{sample_name}.sample', plas_name]) + '\n')
      for sgrna_seq in (sorted(lib.keys())):
        for sgrna_id in lib[sgrna_seq]['ids']:
          count = sample.get(sgrna_id, 0)
          plasmid_count = plasmid.get(sgrna_id, 0)
          f.write('\t'.join([sgrna_id, targeted_genes[sgrna_id], str(count), str(plasmid_count)]) + '\n')
    else:
      f.write('\t'.join(['sgRNA', 'gene', f'{sample_name}.sample']) + '\n')
      for sgrna_seq in (sorted(lib.keys())):
        for sgrna_id in lib[sgrna_seq]['ids']:
          count = sample.get(sgrna_id, 0)
          f.write('\t'.join([sgrna_id, targeted_genes[sgrna_id], str(count)]) + '\n')


def get_singe_guide_library(lib_file: str):
  lib = {}
  targeted_genes = {}

  with open(lib_file, 'r') as f:
    for line_count, line in enumerate(f):
      ele = line.strip().split(',')
      if len(ele) < 3:
        sys.exit(error_msg(f'Guide RNA library file line: {line_count + 1} does not have 3 columns, or the file uses execpted delimiter.'))
      sgrna_id, gene_name, lib_seq = ele[0], ele[1], ele[2]
      if not DNA_PATTERN.match(lib_seq):
        sys.exit(error_msg(f'Sequence column contains non-DNA characters on line: {line_count}.'))

      # NOTE: awkard data structure.
      if lib_seq in lib:
        lib[lib_seq]['ids'].append(sgrna_id)
      else:
        lib[lib_seq] = {'ids': [sgrna_id]}

      targeted_genes[sgrna_id] = gene_name

  return lib, targeted_genes


def get_plasmid_read_counts(plasmid_file: str):
  plasmid, plasmid_name = {}, {}
  with open(plasmid_file, 'r') as f:
    for line_count, line in enumerate(f):
      ele = line.strip().split('\t')
      if len(ele) < 3:
        sys.exit(error_msg(f'Plasmid count file line: {line_count + 1} does not have 3 columns, or the file uses execpted delimiter.'))
      if PLASMID_COUNT_HEADER.match(line):
        plasmid_name = ele[2]
      else:
        sgrna_id, count = ele[0], ele[2]
        plasmid[sgrna_id] = count

  return plasmid, plasmid_name


def get_sgrna_library_counts(input_file: str, ref_file: str, lib: Dict[str, Any], trim: int, reverse_complementing: bool):
  seen = {}

  try:
    samfile = pysam.AlignmentFile(input_file, "rc", reference_filename=ref_file)
  except Exception as e:
    sys.exit(error_msg('Unexpected exception when trying to open input CRAM file: %s' % str(e)))

  sample_name = None
  for rg in samfile.header.to_dict().get('RG'): # does not matter which RG line's SM tag is used
    sample_name = rg.get('SM')
  if not sample_name:
    sys.exit(error_msg('Could not find "SM" tag in the input file header'))

  lib_seqs = {}
  for seq in lib.keys():
    key = seq
    if reverse_complementing:
      # reverse complementing guide RNA sequences instead of each read
      key = rev_compl(key)
    lib_seqs[key] = seq

  # assume library sequences are in same length
  for seq in lib_seqs.keys():
    lib_seq_size = len(seq)
    break

  for read in samfile.fetch(until_eof=True):
    # if the alignment is secondary, supplymentary or vendor failed, skip it!
    if read.flag & 2816:
      continue

    if reverse_complementing:
      cram_seq = read.get_forward_sequence()[-trim-lib_seq_size:-trim]
    else:
      cram_seq = read.get_forward_sequence()[trim:trim+lib_seq_size]

    matching_lib_seq = lib_seqs.get(cram_seq)
    if matching_lib_seq:
      for grna_id in lib[matching_lib_seq]['ids']:
        seen[grna_id] = seen.get(grna_id, 0) + 1

  return seen, sample_name


def count_dual(args):

  #library file must have the following columns defined :
  # unique_id, target_id, gener_pair_id, sgrna_left_seq_id, sgrna_left_seg, sgrna_right_seq_id, sgrna_right_seg
  # unique_id, target_id, gener_pair_id are informative fields that get passed along to output reports

  #Create various lookup dictionaries from the library file

  lookupGuidePair, lookupGuideLeft, lookupGuideRight, lookupGuideLeftRC, lookupGuideRightRC, lookupSafe = {}, {}, {}, {}, {}, {} 
  with open(args['library']) as f:
    header = f.readline()
    fields = header.strip().split('\t')
    #locate the columns in the guide library file
    header_index = {}
    for index, col_name in enumerate(fields):
      for expected_col_name in DUAL_LIBRARY_EXPECTED_HEADER:
        if expected_col_name == col_name.lower():
          header_index[expected_col_name] = index
          break
    # check if all expeted headers are foundi in the input library file
    for expected_col_name in DUAL_LIBRARY_EXPECTED_HEADER:
      if expected_col_name not in header_index.keys():
        sys.exit(error_msg(f'Cound not find named column: {expected_col_name} in the input library file, please check file columns and try again.'))
    
    for line in f:
      eles = line.strip().split('\t')
      sgSeqL = eles[header_index['sgrna_left_seq']]
      sgSeqR = eles[header_index['sgrna_right_seq']]
      sgSeqLrc = rev_compl(sgSeqL)
      sgSeqRrc = rev_compl(sgSeqR)
      lookupGuideLeft[sgSeqL] = 0
      lookupGuideRight[sgSeqR] = 0
      lookupGuideLeftRC[sgSeqLrc]  = 0
      lookupGuideRightRC[sgSeqRrc] = 0
      #store the safe sequences (guide id starts with F followed by a number)
      if SAFE_SEQ_FORMAT.match(eles[header_index['sgrna_left_id']]):
        lookupSafe[sgSeqL] = 0 
      if SAFE_SEQ_FORMAT.match(eles[header_index['sgrna_right_id']]):
        lookupSafe[sgSeqR] = 0 
      lookupGuidePair[sgSeqLrc + sgSeqR] = 0

  n_safe_safe, n_grna1_safe, n_safe_grna2, n_grna1_grna2 = 0, 0, 0, 0
  read_id, n_grna1, n_grna2, n_incorrect_pair, n_miss_miss = None, 0, 0, 0, 0

  with open_plain_or_gzipped_file(
    args['fastq1']) as fq1, open_plain_or_gzipped_file(args['fastq2']) as fq2, open(args['reads'], 'w') as classified_reads:
    for line_index, r1 in enumerate(fq1):
      r2 = fq2.readline()
      residue = (line_index+1)%4  # to figure which of the 4 line of a read recored this line is
      if residue == 1:
        read_id = r1[1:-3]
      elif residue == 2:
        r1 = r1.strip()
        r2 = r2.strip()
        pair_guide = r2 + r1
        # look for correctly paired reads:
        # Reverse Complement (Read2) -> gRNA1 (left); Read1 -> gRNA2 (right)
        if pair_guide in lookupGuidePair:
          lookupGuidePair[pair_guide] += 1
          r2rc = rev_compl(r2)
          label1 = 'safe' if r2rc in lookupSafe else 'gRNA1'
          label2 = 'safe' if r1 in lookupSafe else 'gRNA2'
          #count number of occurrances
          if label1 == 'gRNA1' and label2 == 'gRNA2':
            n_grna1_grna2 += 1
          elif label1 == 'gRNA1' and label2 == 'safe':
            n_grna1_safe += 1
          elif label1 == 'safe' and label2 == 'gRNA2':
            n_safe_grna2 += 1
          else:
            n_safe_safe += 1
          classified_reads.write('\t'.join(['FOUND', f'{label1}_{label2}', args['sample'], read_id, r1, r2, f'{r2rc}{r1}']) + '\n')

        #both guides found but they are incorrectly paired (most reads fall here)
        elif r2 in lookupGuideLeftRC and r1 in lookupGuideRight:
          n_incorrect_pair += 1
          classified_reads.write('\t'.join(['MISS', 'gRNA1_gRNA2', args['sample'], read_id, r1, r2, 'NA']) + '\n')
        # both guides found but they are incorrectly paired and have wrong orientation
        # few reads fall here
        elif r1 in lookupGuideLeft and r2 in lookupGuideRightRC:
          n_incorrect_pair += 1
          classified_reads.write('\t'.join(['MISS', 'gRNA1_gRNA2', args['sample'], read_id, r1, r2, 'NA']) + '\n')
        #only found the left guide (with either correct or wrong orientation)
        elif r2 in lookupGuideLeftRC or r1 in lookupGuideLeft:
          n_grna1 += 1
          classified_reads.write('\t'.join(['MISS', 'gRNA1_nothing', args['sample'], read_id, r1, r2, 'NA']) + '\n')
        #only found the right guide (with either correct or wrong orientation)
        elif r1 in lookupGuideRight or r2 in lookupGuideRightRC:
          n_grna2 += 1
          classified_reads.write('\t'.join(['MISS', 'nothing_gRNA2', args['sample'], read_id, r1, r2, 'NA']) + '\n')
        #didn't match any guides
        else:
          n_miss_miss += 1
          classified_reads.write('\t'.join(['MISS', 'nothing_nothing', args['sample'], read_id, r1, r2, 'NA']) + '\n')

  read_counts = (line_index + 1)/4

  total_guides, zero_guide, less_30_guides = 0, 0, 0
  with open(args['library'], 'r') as lib, open(args['counts'], 'w') as out_count:
    next(lib)
    out_count.write('\t'.join(['unique_id', 'target_id', 'gene_pair_id', 'sample_name']) + '\n')
    for lib_line_index, line in enumerate(lib):
      ele = line.strip().split('\t')
      sgSeqL =ele[header_index['sgrna_left_seq']]
      sgSeqR = ele[header_index['sgrna_right_seq']]
      gene_pair_id = ele[header_index['gene_pair_id']]
      unique_pair_id = ele[header_index['unique_id']]
      target_pair_id = ele[header_index['target_id']]
      sgSeqLrc = rev_compl(sgSeqL)
      counts =lookupGuidePair[sgSeqLrc + sgSeqR]
      out_count.write( '\t'.join([unique_pair_id, target_pair_id, gene_pair_id, str(counts)]) + '\n')
      if counts == 0:
        zero_guides += 1
      if counts < 30:
        less_30_guides += 1

  total_guides = lib_line_index + 1
  with open(args['stats'], 'w') as stats_out:
    stats_out.write(
      '\t'.join(
        ['sample', 'total_reads', 'miss', 'mismatch', 'gRNA1_hits', 'gRNA2_hits',
          'safe_safe', 'gRNA1_safe', 'safe_gRNA2', 'gRNA1_gRNA2', 'total_guides', 'zero_guides', 'less_30_guides']
      ) + '\n')
    stats_out.write('\t'.join(
      [
        args['sample'],
        *[
          str(int(number)) for number in
          [read_counts, n_miss_miss, n_incorrect_pair, n_grna1, n_grna2,
            n_safe_safe, n_grna1_safe, n_safe_grna2, n_grna1_grna2, total_guides, zero_guides, less_30_guides]
        ]
      ]
    ) + '\n')

def rev_compl(dna: str):
  trans_dict = {65: 84, 67: 71, 71: 67, 84: 65, 97: 116, 99: 103, 103: 99, 116: 97}
  return ''.join(reversed(dna)).translate(trans_dict)
