from crispr_read_counts.dual_guide_count import count_dual
import os
import tempfile
import filecmp

test_data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'test-dual')


def test_dual_guide_count():
  args = {
    'library': os.path.join(test_data_dir, 'library_parsed_library_for_counting_without_uveal.test.tsv'),
    'fastq1': os.path.join(test_data_dir, 'A375_c9_day_28_1000x_3_r1.test.fq.gz'),
    'fastq2': os.path.join(test_data_dir, 'A375_c9_day_28_1000x_3_r2.test.fq.gz'),
    'sample': 'test_sample'
  }
  compare_to = {
    'reads': os.path.join(test_data_dir, 'test_dual_classified_reads.test.txt'),
    'stats': os.path.join(test_data_dir, 'test_dual_stats.test.txt'),
    'counts': os.path.join(test_data_dir, 'test_dual_counts.test.txt')
  }
  with tempfile.TemporaryDirectory() as tmpd:
    args['reads'] = os.path.join(tmpd, 'test_dual_classified_reads.test.txt')
    args['stats'] = os.path.join(tmpd, 'test_dual_stats.test.txt')
    args['counts'] = os.path.join(tmpd, 'test_dual_counts.test.txt')
    count_dual(args)
    for option, pointing_file in compare_to.items():
      assert filecmp.cmp(args[option], pointing_file)
