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

  def get_sgrna_library_counts(self, trim: int, reverse_complementing: bool, with_stats: bool):
    if with_stats:
      self._get_sgrna_library_counts_with_stats(trim, reverse_complementing)
    else:
      self._get_sgrna_library_counts_without_stats(trim, reverse_complementing)

  def open_cram_and_get_sample_name(self):
    if not self.ref:
      sys.exit(error_msg(f'Reference file must be provided for reading a CRAM file.'))
    try:
      samfile = pysam.AlignmentFile(self.in_file, "rc", reference_filename=self.ref)
    except Exception as e:
      sys.exit(error_msg('Unexpected exception when trying to open input CRAM file: %s' % str(e)))

    for rg in samfile.header.to_dict().get('RG'): # does not matter which RG line's SM tag is used
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

  def _get_sgrna_library_counts_with_stats(self, trim: int, reverse_complementing: bool):

    total_reads, vendor_failed_reads, mapped_to_guide_reads = 0, 0, 0

    samfile = self.open_cram_and_get_sample_name()
    lib_seqs, lib_seq_size = self.get_lib_seq_dict_and_seq_length(reverse_complementing)

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
      if reverse_complementing:
        # NOTE: it could be as simple as this: cram_seq = read.get_forward_sequence()[-trim-lib_seq_size:-trim]
        # but when trim is 0, end index will be "-0", and this will upset python and return an empty string.
        # so, trim is always added by 1, and the substring became a bit complicated like below.
        cram_seq = read.get_forward_sequence()
        cram_seq = cram_seq[-(trim+1)-(lib_seq_size-1):-(trim+1)] + cram_seq[-(trim+1)]
      else:
        cram_seq = read.get_forward_sequence()[trim:trim+lib_seq_size]

      matching_lib_seq = lib_seqs.get(cram_seq)
      if matching_lib_seq:
        mapped_to_guide_reads += 1
        for grna_id in self.lib[matching_lib_seq]:
          self.sample_count[grna_id] = self.sample_count.get(grna_id, 0) + 1

    self.stats['total_reads'] = total_reads
    self.stats['vendor_failed_reads'] = vendor_failed_reads
    self.stats['mapped_to_guide_reads'] = mapped_to_guide_reads

  def _get_sgrna_library_counts_without_stats(self, trim: int, reverse_complementing: bool):

    samfile = self.open_cram_and_get_sample_name()
    lib_seqs, lib_seq_size = self.get_lib_seq_dict_and_seq_length(reverse_complementing)

    for read in samfile.fetch(until_eof=True):
      # if the alignment is secondary, supplymentary or  vendor failed, skip it!
      if read.flag & 2816:
        continue

      if reverse_complementing:
        # NOTE: it could be as simple as this: cram_seq = read.get_forward_sequence()[-trim-lib_seq_size:-trim]
        # but when trim is 0, end index will be "-0", and this will upset python and return an empty string.
        # so, trim is always added by 1, and the substring became a bit complicated like below.
        cram_seq = read.get_forward_sequence()
        cram_seq = cram_seq[-(trim+1)-(lib_seq_size-1):-(trim+1)] + cram_seq[-(trim+1)]
      else:
        cram_seq = read.get_forward_sequence()[trim:trim+lib_seq_size]

      matching_lib_seq = lib_seqs.get(cram_seq)
      if matching_lib_seq:
        for grna_id in self.lib[matching_lib_seq]:
          self.sample_count[grna_id] = self.sample_count.get(grna_id, 0) + 1

  def write_output(self, out_stats: str):
    if out_stats:
      self._write_output_with_stats(out_stats)
    else:
      self._write_output_without_stats()

  def _write_output_with_stats(self, out_stats: str):
    zero_count_guides, low_count_guides = 0, 0
    with open(self.out_count, 'w') as f:
      if self.plas_name:
        f.write('\t'.join(['sgRNA', 'gene', f'{self.sample_name}.sample', self.plas_name]) + '\n')
        for sgrna_seq in (sorted(self.lib.keys())):
          for sgrna_id in self.lib[sgrna_seq]:
            count = self.sample_count.get(sgrna_id, 0)
            if count == 0:
              zero_count_guides += 1
            if count < self.LOW_COUNT_GUIDES_THRESHOLD:
              low_count_guides += 1
            plasmid_count = self.plasmid.get(sgrna_id, 0)
            f.write('\t'.join([sgrna_id, self.targeted_genes[sgrna_id], str(count), str(plasmid_count)]) + '\n')
      else:
        f.write('\t'.join(['sgRNA', 'gene', f'{self.sample_name}.sample']) + '\n')
        for sgrna_seq in (sorted(self.lib.keys())):
          for sgrna_id in self.lib[sgrna_seq]:
            count = self.sample_count.get(sgrna_id, 0)
            if count == 0:
              zero_count_guides += 1
            if count < self.LOW_COUNT_GUIDES_THRESHOLD:
              low_count_guides += 1
            f.write('\t'.join([sgrna_id, self.targeted_genes[sgrna_id], str(count)]) + '\n')

    self.stats['zero_count_guides'] = zero_count_guides
    self.stats['low_count_guides'] = low_count_guides
    with open(out_stats, 'w') as out_s:
      json.dump(self.stats, out_s)
      out_s.write('\n')

  def _write_output_without_stats(self):
    with open(self.out_count, 'w') as f:
      if self.plas_name:
        f.write('\t'.join(['sgRNA', 'gene', f'{self.sample_name}.sample', self.plas_name]) + '\n')
        for sgrna_seq in (sorted(self.lib.keys())):
          for sgrna_id in self.lib[sgrna_seq]:
            count = self.sample_count.get(sgrna_id, 0)
            plasmid_count = self.plasmid.get(sgrna_id, 0)
            f.write('\t'.join([sgrna_id, self.targeted_genes[sgrna_id], str(count), str(plasmid_count)]) + '\n')
      else:
        f.write('\t'.join(['sgRNA', 'gene', f'{self.sample_name}.sample']) + '\n')
        for sgrna_seq in (sorted(self.lib.keys())):
          for sgrna_id in self.lib[sgrna_seq]:
            count = self.sample_count.get(sgrna_id, 0)
            f.write('\t'.join([sgrna_id, self.targeted_genes[sgrna_id], str(count)]) + '\n')

  def count(self, trim, plasmid_count_file, reverse_complement, out_stats):
    if plasmid_count_file:
      self.plasmid, self.plas_name = self.get_plasmid_read_counts(plasmid_count_file)
    self.get_sgrna_library_counts(trim, reverse_complement, bool(out_stats))
    self.write_output(out_stats)

  @staticmethod
  def get_single_guide_library(lib_file: str):
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

        if lib_seq in lib:
          lib[lib_seq].append(sgrna_id)
        else:
          lib[lib_seq] = [sgrna_id]

        targeted_genes[sgrna_id] = gene_name

    return lib, targeted_genes

  @staticmethod
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
