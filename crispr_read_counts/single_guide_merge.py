import sys
import json
from .utils import (
  error_msg,
  open_plain_or_gzipped_file,
  check_file_readable,
  check_file_writable,
  PLASMID_COUNT_HEADER)
from typing import TextIO, List, Dict
from .single_guide_count import SingleGuideReadCounts


def merge_single(args):
  '''
  handler of user inputs and execute the workflow
  '''
  has_plasmid = args['plasmid']
  files = args['input'].split(',')
  for an_input in files:
    check_file_readable(an_input, f'Provided read counts file does not exist or have no permission to read: {an_input}')
  if args['stats']:
    check_file_writable(args['stats'], 'Cannot write to provided output stats file: %s' % args['stats'])

  samp_name, plas_name, sample_rc, plasmid_rc, genes = get_sample_read_counts(files, has_plasmid)
  print(f'writing merged counts to: {args["output"]}...', flush=True)
  with open(args['output'], 'w', newline='') as out:
    if has_plasmid:
      out.write('\t'.join(['sgRNA', 'gene', samp_name, plas_name]) + '\n')
      for id in sample_rc.keys():
        out.write('\t'.join([id, genes[id], str(sample_rc[id]), str(plasmid_rc[id])]) + '\n')
    else:
      out.write('\t'.join(['sgRNA', 'gene', samp_name]) + '\n')
      for id in sample_rc.keys():
        out.write('\t'.join([id, genes[id], str(sample_rc[id])]) + '\n')

  if args['stats']:
    print(f'writing stats to: {args["stats"]}...', flush=True)
    write_stats_to_file(sample_rc, args['stats'])
  print('Done.')


def get_sample_read_counts(files: List[str], has_plasmid: bool):
  sample_name, plasmid_name = None, None
  sample, plasmid, targeted_genes = {}, {}, {}

  for a_file in files:
    print(f'reading from {a_file}...')
    with open_plain_or_gzipped_file(a_file) as in_f:
      header = in_f.readline().strip()
      header_split = header.split('\t')
      if PLASMID_COUNT_HEADER.match(header):
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

      def get_counts_results(in_f: TextIO, with_or_without_plasmid_counts: bool):
        '''
        read a count file and update sample dict, targeted_genes dict and plasmid dict if plasmid counts are available.
        '''
        if with_or_without_plasmid_counts:
          for line in in_f:
            line_split = line.strip().split('\t')
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
          for line in in_f:
            line_split = line.strip().split('\t')
            id = line_split[0]
            sample_count = sample.get(id, 0) + int(line_split[2])
            sample[id] = sample_count
            targeted_genes[id] = line_split[1]

      get_counts_results(in_f, has_plasmid)

  return sample_name, plasmid_name, sample, plasmid, targeted_genes


def write_stats_to_file(sample_count: Dict[str, int], output_path: str):
  stats = {
    'zero_count_guides': 0,
    'low_count_guides': 0,
    'total_counts': 0
  }

  for a_count in sample_count.values():
    stats['total_counts'] += a_count
    if a_count == 0:
      stats['zero_count_guides'] += 1
    if a_count < SingleGuideReadCounts.LOW_COUNT_GUIDES_THRESHOLD:
      stats['low_count_guides'] += 1
  with open(output_path, 'w') as out_s:
    json.dump(stats, out_s)
    out_s.write('\n')
