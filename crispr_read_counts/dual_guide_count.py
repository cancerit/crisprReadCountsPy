
import sys
from typing import List
from .utils import (
  error_msg,
  warning_msg,
  open_plain_or_gzipped_file,
  rev_compl,
  SAFE_SEQ_FORMAT,
  check_file_readable,
  check_file_writable)

DUAL_LIBRARY_EXPECTED_HEADER = ['sgrna_left_id', 'sgrna_left_seq', 'sgrna_right_id', 'sgrna_right_seq', 'unique_id', 'gene_pair_id', 'target_id']


def count_dual(args):

  # library file must have the following columns defined :
  # unique_id, target_id, gener_pair_id, sgrna_left_seq_id, sgrna_left_seg, sgrna_right_seq_id, sgrna_right_seg
  # unique_id, target_id, gener_pair_id are informative fields that get passed along to output reports

  # Create various lookup dictionaries from the library file
  validate_inputs(args)
  (lookupGuidePair, lookupGuideLeft, lookupGuideRight, lookupGuideLeftRC,
   lookupGuideRightRC, lookupSafe, header_index) = library_to_dicts(args['library'])

  (n_safe_safe, n_grna1_safe, n_safe_grna2,
   n_grna1_grna2, n_grna1, n_grna2, n_incorrect_pair, n_miss_miss, read_counts
   ) = write_classified_reads_to_file_return_stats(
      args['fastq1'], args['fastq2'], args['reads'], args['sample'],
      lookupGuidePair, lookupGuideLeft, lookupGuideRight, lookupGuideLeftRC, lookupGuideRightRC, lookupSafe)

  total_guides, zero_guides, less_30_guides = write_guides_return_stats(
    args['library'], args['counts'], args['sample'], lookupGuidePair, header_index)

  write_stats(
    args['stats'],
    ['sample', 'total_reads', 'miss', 'mismatch', 'gRNA1_hits', 'gRNA2_hits', 'safe_safe',
     'gRNA1_safe', 'safe_gRNA2', 'gRNA1_gRNA2', 'total_guides', 'zero_guides', 'less_30_guides'],
    [
      args['sample'],
      *[
        str(int(number)) for number in
        [read_counts, n_miss_miss, n_incorrect_pair, n_grna1, n_grna2,
         n_safe_safe, n_grna1_safe, n_safe_grna2, n_grna1_grna2, total_guides, zero_guides, less_30_guides]
      ]
    ])


def validate_inputs(args):
  for file_type, file_path in zip(['library', 'FastQ', 'FastQ'], [args['library'], args['fastq1'], args['fastq2']]):
    check_file_readable(file_path, f'Provided {file_type} file does not exist or have no permission to read: {file_path}')

  for file_type, file_path in zip(['classified reads', 'counts', 'stats'], [args['reads'], args['counts'], args['stats']]):
    check_file_writable(file_path, f'Cannot write to provided output {file_type} file: {file_path}.')


def library_to_dicts(library: str):

  lookupGuidePair, lookupGuideLeft, lookupGuideRight, lookupGuideLeftRC, lookupGuideRightRC, lookupSafe, header_index = {}, {}, {}, {}, {}, {}, {}

  with open(library) as f:
    header = f.readline().strip().split('\t')
    # locate the columns in the guide library file
    for index, col_name in enumerate(header):
      for expected_col_name in DUAL_LIBRARY_EXPECTED_HEADER:
        if expected_col_name == col_name.lower():
          header_index[expected_col_name] = index
          break
    # check if all expeted headers are foundi in the input library file
    for expected_col_name in DUAL_LIBRARY_EXPECTED_HEADER:
      if expected_col_name not in header_index.keys():
        sys.exit(error_msg(f'Cound not find named column: {expected_col_name} in the input library file, please check file columns and try again.'))

    for line in f:
      line_split = line.strip().split('\t')
      sgSeqL = line_split[header_index['sgrna_left_seq']]
      sgSeqR = line_split[header_index['sgrna_right_seq']]
      sgSeqLrc = rev_compl(sgSeqL)
      sgSeqRrc = rev_compl(sgSeqR)
      lookupGuideLeft[sgSeqL] = 0
      lookupGuideRight[sgSeqR] = 0
      lookupGuideLeftRC[sgSeqLrc] = 0
      lookupGuideRightRC[sgSeqRrc] = 0
      # store the safe sequences (guide id starts with F followed by a number)
      if SAFE_SEQ_FORMAT.match(line_split[header_index['sgrna_left_id']]):
        lookupSafe[sgSeqL] = 0
      if SAFE_SEQ_FORMAT.match(line_split[header_index['sgrna_right_id']]):
        lookupSafe[sgSeqR] = 0
      lookupGuidePair[sgSeqLrc + sgSeqR] = 0

  return lookupGuidePair, lookupGuideLeft, lookupGuideRight, lookupGuideLeftRC, lookupGuideRightRC, lookupSafe, header_index


def write_classified_reads_to_file_return_stats(
  fastq1: str, fastq2: str, out_reads: str, sample_name: str,
  lookupGuidePair, lookupGuideLeft, lookupGuideRight, lookupGuideLeftRC, lookupGuideRightRC, lookupSafe):

  n_safe_safe, n_grna1_safe, n_safe_grna2, n_grna1_grna2, n_grna1, n_grna2, n_incorrect_pair, n_miss_miss = 0, 0, 0, 0, 0, 0, 0, 0

  read_id = None
  with open_plain_or_gzipped_file(
    fastq1) as fq1, open_plain_or_gzipped_file(fastq2) as fq2, open(out_reads, 'w') as classified_reads:
    for line_index, r1 in enumerate(fq1, 1):
      r2 = fq2.readline()
      residue = (line_index) % 4  # to figure which of the 4 line of a read recored this line is
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
          # count number of occurrances
          if label1 == 'gRNA1' and label2 == 'gRNA2':
            n_grna1_grna2 += 1
          elif label1 == 'gRNA1' and label2 == 'safe':
            n_grna1_safe += 1
          elif label1 == 'safe' and label2 == 'gRNA2':
            n_safe_grna2 += 1
          else:
            n_safe_safe += 1
          classified_reads.write('\t'.join(['FOUND', f'{label1}_{label2}', sample_name, read_id, r1, r2, f'{r2rc}{r1}']) + '\n')

        # both guides found but they are incorrectly paired (most reads fall here)
        elif r2 in lookupGuideLeftRC and r1 in lookupGuideRight:
          n_incorrect_pair += 1
          classified_reads.write('\t'.join(['MISS', 'gRNA1_gRNA2', sample_name, read_id, r1, r2, 'NA']) + '\n')
        # both guides found but they are incorrectly paired and have wrong orientation
        # few reads fall here
        elif r1 in lookupGuideLeft and r2 in lookupGuideRightRC:
          n_incorrect_pair += 1
          classified_reads.write('\t'.join(['MISS', 'gRNA1_gRNA2', sample_name, read_id, r1, r2, 'NA']) + '\n')
        # only found the left guide (with either correct or wrong orientation)
        elif r2 in lookupGuideLeftRC or r1 in lookupGuideLeft:
          n_grna1 += 1
          classified_reads.write('\t'.join(['MISS', 'gRNA1_nothing', sample_name, read_id, r1, r2, 'NA']) + '\n')
        # only found the right guide (with either correct or wrong orientation)
        elif r1 in lookupGuideRight or r2 in lookupGuideRightRC:
          n_grna2 += 1
          classified_reads.write('\t'.join(['MISS', 'nothing_gRNA2', sample_name, read_id, r1, r2, 'NA']) + '\n')
        # didn't match any guides
        else:
          n_miss_miss += 1
          classified_reads.write('\t'.join(['MISS', 'nothing_nothing', sample_name, read_id, r1, r2, 'NA']) + '\n')

  if (line_index) % 4 != 0:
    print(warning_msg('Number of lines in provided FastQ files is not multiple times of 4, truncated file?'), flush=True)

  read_counts = int((line_index + 1) / 4)

  return n_safe_safe, n_grna1_safe, n_safe_grna2, n_grna1_grna2, n_grna1, n_grna2, n_incorrect_pair, n_miss_miss, read_counts


def write_guides_return_stats(library: str, out_counts: str, sample_name: str, lookupGuidePair, header_index):
  zero_guides, less_30_guides = 0, 0
  with open(library, 'r') as lib, open(out_counts, 'w') as out_ct:
    next(lib)
    out_ct.write('\t'.join(['unique_id', 'target_id', 'gene_pair_id', sample_name]) + '\n')
    for line_count, line in enumerate(lib, 1):
      ele = line.strip().split('\t')
      sgSeqL = ele[header_index['sgrna_left_seq']]
      sgSeqR = ele[header_index['sgrna_right_seq']]
      gene_pair_id = ele[header_index['gene_pair_id']]
      unique_pair_id = ele[header_index['unique_id']]
      target_pair_id = ele[header_index['target_id']]
      sgSeqLrc = rev_compl(sgSeqL)
      counts = lookupGuidePair[sgSeqLrc + sgSeqR]
      out_ct.write('\t'.join([unique_pair_id, target_pair_id, gene_pair_id, str(counts)]) + '\n')
      if counts == 0:
        zero_guides += 1
      if counts < 30:
        less_30_guides += 1

  total_guides = line_count

  return total_guides, zero_guides, less_30_guides


def write_stats(out_stats: str, col_names: List[str], values: List[str]):

  with open(out_stats, 'w', newline='') as stats_out:
    stats_out.write('\t'.join(col_names) + '\n')
    stats_out.write('\t'.join(values) + '\n')
