import sys
from typing import Dict, Any
from .utils import (
  error_msg,
  open_plain_or_gzipped_file,
  rev_compl,
  PLASMID_COUNT_HEADER,
  DNA_PATTERN,
  check_file_readable,
  check_file_writable)
import pysam
import json
import csv

csv.register_dialect('sg_library', delimiter=',', quoting=csv.QUOTE_NONE)
csv.register_dialect('sg_plasmid_count', delimiter='\t', quoting=csv.QUOTE_NONE)
csv.register_dialect('sg_out_count', delimiter='\t', quoting=csv.QUOTE_NONE)


def count_single(args: Dict[str, Any]):
  # validate inputs before doing anything
  check_input_files(args)
  count_instance = SingleGuideReadCounts(args['library'], args['input'], args['output'], args['ref'])
  count_instance.count(args['trim'], args['plasmid'], args['reverse_complement'], args['stats'])


def check_input_files(args: Dict[str, Any]):
  if not check_file_readable(args['library']):
    sys.exit(error_msg('Provided library file does not exist or have no permission to read: %s' % {args['library']}))
  if not check_file_readable(args['input']):
    sys.exit(error_msg('Provided input file does not exist or have no permission to read: %s' % args['input']))
  if args['ref'] and not check_file_readable(args['ref']):
    sys.exit(error_msg('Provided reference file does not exist or have no permission to read: %s' % args['ref']))
  if args['plasmid'] and not check_file_readable(args['plasmid']):
    sys.exit(error_msg('Provided plasmid count file does not exist or have no permission to read: %s' % args['plasmid']))
  if not check_file_writable(args['output']):
    sys.exit(error_msg('Cannot write to provided output count file: %s' % args['output']))
  if args['stats'] and not check_file_writable(args['stats']):
    sys.exit(error_msg('Cannot write to provided output stats file: %s' % args['stats']))


class SingleGuideReadCounts():
  '''
  The class is just to reduce parameters passing around functions.
  '''

  LOW_COUNT_GUIDES_THRESHOLD = 15

  def __init__(self, library, in_file, out_count, ref):
    self.lib, self.targeted_genes = self.get_single_guide_library(library)
    self.in_file = in_file
    self.out_count = out_count
    self.ref = ref  # ref must be an existing file if input is a CRAM
    self.plasmid = None
    self.plas_name = None
    self.sample_count = {}
    self.sample_name = None
    self.stats = {}

  def open_cram_and_get_sample_name(self):
    if not self.ref:
      sys.exit(error_msg(f'Reference file must be provided for reading a CRAM file.'))
    try:
      samfile = pysam.AlignmentFile(self.in_file, "rc", reference_filename=self.ref)
    except Exception as e:
      sys.exit(error_msg('Unexpected exception when trying to open input CRAM file: %s' % str(e)))

    for rg in samfile.header.to_dict().get('RG'):  # does not matter which RG line's SM tag is used
      self.sample_name = rg.get('SM')
    if not self.sample_name:
      sys.exit(error_msg('Could not find "SM" tag in the input file header'))
    return samfile

  def get_lib_seq_dict_and_seq_length(self, reverse_complementing):
    lib_seqs = {}
    for seq in self.lib.keys():
      key = seq
      if reverse_complementing:
        # reverse complementing guide RNA sequences instead of each read
        key = rev_compl(key)
      lib_seqs[key] = seq

    # assume library sequences are in same length
    for seq in lib_seqs.keys():
      lib_seq_size = len(seq)
      break

    return lib_seqs, lib_seq_size

  def get_sgrna_library_counts(self, trim: int, reverse_complementing: bool):
    '''
    # NOTE: Stats are calculated regardless whether they're required or not in order to achieve better code maintainability.
    # From limited benchmarking runs, this only increase ~2% run time with 11 million reads as input.
    '''
    total_reads, vendor_failed_reads, mapped_to_guide_reads = 0, 0, 0

    samfile = self.open_cram_and_get_sample_name()
    lib_seqs, lib_seq_size = self.get_lib_seq_dict_and_seq_length(reverse_complementing)
    sl = self.get_seq_slicing_indexes(reverse_complementing, trim, lib_seq_size)

    for read in samfile.fetch(until_eof=True):
      # if the alignment is secondary or supplymentary, skip it!
      if read.flag & 2304:
        continue

      # if the alignment is vendor failed, skip it but count it!
      if read.flag & 512:
        total_reads += 1
        vendor_failed_reads += 1
        continue

      total_reads += 1
      cram_seq = read.get_forward_sequence()[sl]

      matching_lib_seq = lib_seqs.get(cram_seq)
      if matching_lib_seq:
        mapped_to_guide_reads += 1
        for grna_id in self.lib[matching_lib_seq]:
          self.sample_count[grna_id] = self.sample_count.get(grna_id, 0) + 1

    self.stats['total_reads'] = total_reads
    self.stats['vendor_failed_reads'] = vendor_failed_reads
    self.stats['mapped_to_guide_reads'] = mapped_to_guide_reads

  def write_output(self, out_stats: str):
    zero_count_guides, low_count_guides = 0, 0
    with open(self.out_count, 'w', newline='') as f:
      writer = csv.writer(f, 'sg_out_count')
      if self.plas_name:
        writer.writerow(['sgRNA', 'gene', f'{self.sample_name}.sample', self.plas_name])
        for sgrna_seq in sorted(self.lib.keys()):
          for sgrna_id in self.lib[sgrna_seq]:
            count = self.sample_count.get(sgrna_id, 0)
            if count == 0:
              zero_count_guides += 1
            if count < self.LOW_COUNT_GUIDES_THRESHOLD:
              low_count_guides += 1
            plasmid_count = self.plasmid.get(sgrna_id, 0)
            writer.writerow([sgrna_id, self.targeted_genes[sgrna_id], str(count), str(plasmid_count)])
      else:
        writer.writerow(['sgRNA', 'gene', f'{self.sample_name}.sample'])
        for sgrna_seq in sorted(self.lib.keys()):
          for sgrna_id in self.lib[sgrna_seq]:
            count = self.sample_count.get(sgrna_id, 0)
            if count == 0:
              zero_count_guides += 1
            if count < self.LOW_COUNT_GUIDES_THRESHOLD:
              low_count_guides += 1
            writer.writerow([sgrna_id, self.targeted_genes[sgrna_id], str(count)])

    if out_stats:
      self.stats['zero_count_guides'] = zero_count_guides
      self.stats['low_count_guides'] = low_count_guides
      with open(out_stats, 'w') as out_s:
        json.dump(self.stats, out_s)
        out_s.write('\n')

  def count(self, trim, plasmid_count_file, reverse_complement, out_stats):
    if plasmid_count_file:
      self.plasmid, self.plas_name = self.get_plasmid_read_counts(plasmid_count_file)
    self.get_sgrna_library_counts(trim, reverse_complement)
    self.write_output(out_stats)

  @staticmethod
  def get_single_guide_library(lib_file: str):
    lib = {}
    targeted_genes = {}

    with open(lib_file, 'r') as f:
      reader = csv.reader(f, 'sg_library')
      for line_number, line_split in enumerate(reader, 1):
        if len(line_split) < 3:
          sys.exit(error_msg(f'Guide RNA library file line: {line_number} does not have 3 columns, or the file uses expected delimiter.'))
        sgrna_id, gene_name, lib_seq = line_split[0], line_split[1], line_split[2]
        if not DNA_PATTERN.match(lib_seq):
          sys.exit(error_msg(f'Sequence column contains non-DNA characters on line: {line_number}.'))

        if lib_seq in lib:
          lib[lib_seq].append(sgrna_id)
        else:
          lib[lib_seq] = [sgrna_id]

        targeted_genes[sgrna_id] = gene_name

    return lib, targeted_genes

  @staticmethod
  def get_plasmid_read_counts(plasmid_file: str):
    plasmid = {}
    plasmid_name = None
    with open(plasmid_file, 'r') as f:
      reader = csv.reader(f, 'sg_plasmid_count')
      for line_number, line_split in enumerate(reader, 1):
        if len(line_split) < 3:
          sys.exit(error_msg(f'Plasmid count file line: {line_number} does not have 3 columns, or the file uses expected delimiter.'))
        if line_number == 1:
          if PLASMID_COUNT_HEADER.match('\t'.join(line_split)):
            plasmid_name = line_split[2]
          else:
            sys.exit(error_msg(f'Plasmid count file does not have expected header.'))
        else:
          sgrna_id, count = line_split[0], line_split[2]
          plasmid[sgrna_id] = count

    return plasmid, plasmid_name

  @staticmethod
  def get_seq_slicing_indexes(reverse_complementing, trim, lib_seq_size):
    return (
      slice(trim, trim + lib_seq_size) if not reverse_complementing else
      slice(-trim - lib_seq_size, -trim) if trim > 0 else
      slice(-lib_seq_size, None)
    )
