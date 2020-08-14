import sys
from .utils import error_msg, open_plain_or_gzipped_file, PLASMID_COUNT_HEADER
from typing import TextIO
import csv

csv.register_dialect('sg_out_count', delimiter='\t', quoting=csv.QUOTE_NONE)


def merge_single(args):
  has_plasmid = args['plasmid']
  samp_name, plas_name, sample_rc, plasmid_rc, genes = get_sample_read_counts(args['input'], has_plasmid)
  with open(args['output'], 'w', newline='') as out:
    writer = csv.writer(out, 'sg_out_count')
    if has_plasmid:
      writer.writerow(['sgRNA', 'gene', samp_name, plas_name])
      for id in sample_rc.keys():
        writer.writerow([id, genes[id], str(sample_rc[id]), str(plasmid_rc[id])])
    else:
      writer.writerow(['sgRNA', 'gene', samp_name])
      for id in sample_rc.keys():
        writer.writerow([id, genes[id], str(sample_rc[id])])
  print('Done.')


def get_sample_read_counts(in_files_string: str, has_plasmid: bool):
  sample_name, plasmid_name = None, None
  sample, plasmid, targeted_genes = {}, {}, {}

  files = in_files_string.split(',')
  for a_file in files:
    print(f'reading from {a_file}...')
    with open_plain_or_gzipped_file(a_file) as in_f:
      reader = csv.reader(in_f, 'sg_out_count')
      header_split = reader.__next__()
      if PLASMID_COUNT_HEADER.match('\t'.join(header_split)):
        sample_name = header_split[2]
        if has_plasmid:
          if len(header_split) < 4:
            sys.exit(error_msg(f'Can not find plasmid count column in input file: {a_file}.\nProbably should remove option "--plasmid"?'))
          if plasmid_name is None:
            plasmid_name = header_split[3]
          elif plasmid_name != header_split[3]:
            # files should have same plasmid sample name
            sys.exit(error_msg(f'Plasmid sample names is different in this file: {a_file} from in file: {files[0]}'))
      else:
        sys.exit(error_msg(f'Unexpected header in input file: {a_file}'))

      def get_counts_results(reader: csv.reader, with_or_without_plasmid_counts: bool):
        if with_or_without_plasmid_counts:
          for line_split in reader:
            id = line_split[0]
            sample_count = sample.get(id, 0) + int(line_split[2])
            sample[id] = sample_count
            targeted_genes[id] = line_split[1]
            count = int(line_split[3])
            if id in plasmid and count != plasmid[id]:
              sys.exit(error_msg(f'Plasmid count of sgRNA: {id} is not consistent across input count files.'))
            else:
              plasmid[id] = count
        else:
          for line_split in reader:
            id = line_split[0]
            sample_count = sample.get(id, 0) + int(line_split[2])
            sample[id] = sample_count
            targeted_genes[id] = line_split[1]

      get_counts_results(reader, has_plasmid)

  return sample_name, plasmid_name, sample, plasmid, targeted_genes
