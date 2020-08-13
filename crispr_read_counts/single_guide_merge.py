import sys
from .utils import error_msg, open_plain_or_gzipped_file, PLASMID_COUNT_HEADER
from typing import TextIO


def merge_single(args):
  has_plasmid = args['plasmid']
  samp_name, plas_name, sample_rc, plasmid_rc, genes = get_sample_read_counts(args['input'], has_plasmid)
  with open(args['output'], 'w') as out:
    if has_plasmid:
      out.write('\t'.join(['sgRNA', 'gene', samp_name, plas_name]) + '\n')
      for id in sample_rc.keys():
        out.write('\t'.join([id, genes[id], str(sample_rc[id]), str(plasmid_rc[id])]) + '\n')
    else:
      out.write('\t'.join(['sgRNA', 'gene', samp_name]) + '\n')
      for id in sample_rc.keys():
        out.write('\t'.join([id, genes[id], str(sample_rc[id])]) + '\n')
  print('Done.')


def get_sample_read_counts(in_files_string: str, has_plasmid: bool):
  sample_name, plasmid_name = None, None
  sample, plasmid, targeted_genes = {}, {}, {}

  files = in_files_string.split(',')
  for a_file in files:
    print(f'reading from {a_file}...')
    with open_plain_or_gzipped_file(a_file) as in_f:
      header = in_f.readline()
      if PLASMID_COUNT_HEADER.match(header):
        data = header.strip().split('\t')
        sample_name = data[2]
        if has_plasmid:
          if len(data) < 4:
            sys.exit(error_msg(f'Can not find plasmid count column in input file: {a_file}.\nProbably should remove option "--plasmid"?'))
          if plasmid_name is None:
            plasmid_name = data[3]
          elif plasmid_name != data[3]:
            # files should have same plasmid sample name
            sys.exit(error_msg(f'Plasmid sample names is different in this file: {a_file} from in file: {files[0]}'))
      else:
        sys.exit(error_msg(f'Unexpected header in input file: {a_file}'))

      def get_counts_results(in_f: TextIO, with_or_without_plasmid_counts: bool):
        if with_or_without_plasmid_counts:
          for line in in_f:
            data = line.strip().split('\t')
            id = data[0]
            sample_count = sample.get(id, 0) + int(data[2])
            sample[id] = sample_count
            targeted_genes[id] = data[1]
            if id in plasmid and int(data[3]) != plasmid[id]:
              sys.exit(error_msg(f'Plasmid count of sgRNA: {id} is not consistent across input count files.'))
            else:
              plasmid[id] = int(data[3])
        else:
          for line in in_f:
            data = line.strip().split('\t')
            id = data[0]
            sample_count = sample.get(id, 0) + int(data[2])
            sample[id] = sample_count
            targeted_genes[id] = data[1]

      get_counts_results(in_f, has_plasmid)

  return sample_name, plasmid_name, sample, plasmid, targeted_genes
