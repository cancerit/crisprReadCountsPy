import pytest
from typing import List
from crispr_read_counts.single_guide_count import check_input_files, count_single
from crispr_read_counts.single_guide_merge import merge_single
import os
import tempfile
import filecmp

test_data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
test_single_data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'test-single')

TEST_INPUTS = {'library': os.path.join(test_single_data_dir, 'Human_v1_CRISPR_library.test.lib.csv'),
  'input': os.path.join(test_single_data_dir, 'test.crispr.cram'),
  'output': 'output.txt',
  'stats': 'stats.txt',
  'plasmid': os.path.join(test_single_data_dir, 'Human_v1_CRISPR_library.test.plasmid_counts.csv'),
  'ref': os.path.join(test_single_data_dir, 'genome.empty.fa'),
  'trim': 0,
  'reverse_complement': False,
  'lib_delimiter': ','}


@pytest.mark.parametrize('args', [
  ({'library': os.path.join(test_data_dir, 'non_existing_file')}),
  ({'library': os.path.join(test_single_data_dir, 'Human_v1_CRISPR_library.test.lib.csv'),
    'input': os.path.join(test_single_data_dir, 'non_existing_file')}),
  ({'library': os.path.join(test_single_data_dir, 'Human_v1_CRISPR_library.test.lib.csv'),
    'input': os.path.join(test_single_data_dir, 'test.crispr.cram'),
    'ref': None, 'plasmid': None,
    'output': os.path.join(test_data_dir, 'non_existing_folder', 'non_existing_file')}),
  ({'library': os.path.join(test_single_data_dir, 'Human_v1_CRISPR_library.test.lib.csv'),
    'input': os.path.join(test_single_data_dir, 'test.crispr.cram'),
    'ref': None, 'plasmid': None,
    'output': os.path.join(test_single_data_dir, 'file_to_write_to.txt'),
    'stats': os.path.join(test_data_dir, 'non_existing_folder', 'file_to_write_to.txt')}),
  ({'library': os.path.join(test_single_data_dir, 'Human_v1_CRISPR_library.test.lib.csv'),
    'input': os.path.join(test_single_data_dir, 'test.crispr.cram'),
    'ref': None,
    'output': os.path.join(test_single_data_dir, 'file_to_write_to.txt'),
    'stats': os.path.join(test_single_data_dir, 'file_to_write_to.txt'),
    'plasmid': os.path.join(test_data_dir, 'non_existing_file')}),
  ({'library': os.path.join(test_single_data_dir, 'Human_v1_CRISPR_library.test.lib.csv'),
    'input': os.path.join(test_single_data_dir, 'non_exist_file'),
    'output': os.path.join(test_single_data_dir, 'file_to_write_to.txt'),
    'stats': os.path.join(test_single_data_dir, 'file_to_write_to.txt'),
    'plasmid': os.path.join(test_single_data_dir, 'Human_v1_CRISPR_library.test.plasmid_counts.csv'),
    'ref': os.path.join(test_data_dir, 'non_existing_file')})
])
def test_check_input_files(args: List):
  with pytest.raises(SystemExit) as pytest_e:
    check_input_files(args)
  assert pytest_e.type == SystemExit


@pytest.mark.parametrize('args, expected_output', [
  ({**TEST_INPUTS, 'plasmid': None},
  {'output': os.path.join(test_single_data_dir, 'test.crispr.count.no_plasmid.txt')}),
  ({**TEST_INPUTS, 'plasmid': None,
   'library': os.path.join(test_single_data_dir, 'Human_v1_CRISPR_library.test.lib.tsv'), 'lib_delimiter': '\t'},
  {'output': os.path.join(test_single_data_dir, 'test.crispr.count.no_plasmid.txt')}),
  ({**TEST_INPUTS},
  {'output': os.path.join(test_single_data_dir, 'test.crispr.count.with_plasmid.txt'),
   'stats': os.path.join(test_single_data_dir, 'test.crispr.count.with_plasmid.stats.txt')}),
  ({**TEST_INPUTS, 'trim': 2},
  {'output': os.path.join(test_single_data_dir, 'test.crispr.count.with_plasmid.trim2.txt'),
   'stats': os.path.join(test_single_data_dir, 'test.crispr.count.with_plasmid.trim2.stats.txt')}),
  ({**TEST_INPUTS, 'reverse_complement': True},
  {'output': os.path.join(test_single_data_dir, 'test.crispr.count.with_plasmid.reverse_comp.txt'),
   'stats': os.path.join(test_single_data_dir, 'test.crispr.count.with_plasmid.reverse_comp.stats.txt')}),
  ({**TEST_INPUTS, 'reverse_complement': True, 'trim': 2},
  {'output': os.path.join(test_single_data_dir, 'test.crispr.count.with_plasmid.reverse_comp.trim2.txt'),
   'stats': os.path.join(test_single_data_dir, 'test.crispr.count.with_plasmid.reverse_comp.trim2.stats.txt')}),
  ({**TEST_INPUTS, 'reverse_complement': True, 'trim': 2,
    'library': os.path.join(test_single_data_dir, 'Human_v1_CRISPR_library.test.lib.tsv'), 'lib_delimiter': '\t'},
  {'output': os.path.join(test_single_data_dir, 'test.crispr.count.with_plasmid.reverse_comp.trim2.txt'),
   'stats': os.path.join(test_single_data_dir, 'test.crispr.count.with_plasmid.reverse_comp.trim2.stats.txt')})
])
def test_single_count(args, expected_output):
  with tempfile.TemporaryDirectory() as tmpd:
    args['output'] = os.path.join(tmpd, args['output'])
    args['stats'] = os.path.join(tmpd, args['stats'])
    print(args)
    count_single(args)
    for key, pointing_file in expected_output.items():
      assert filecmp.cmp(args[key], pointing_file)


@pytest.mark.parametrize('args, compare_to', [
  ({'input': '{0},{0}'.format(os.path.join(test_single_data_dir, 'test.crispr.count.with_plasmid.txt')),
    'plasmid': True},
  os.path.join(test_single_data_dir, 'test.crispr.count.with_plasmid.merge_doubled.txt')),
  ({'input': '{0},{0}'.format(os.path.join(test_single_data_dir, 'test.crispr.count.no_plasmid.txt')),
    'plasmid': False},
  os.path.join(test_single_data_dir, 'test.crispr.count.no_plasmid.merge_doubled.txt'))
])
def test_merge_single(args, compare_to):
  with tempfile.TemporaryDirectory() as tmpd:
    args['output'] = os.path.join(tmpd, 'merge_output.txt')
    merge_single(args)
    assert filecmp.cmp(args['output'], compare_to)
