import pytest
from typing import List
from click.testing import CliRunner
from crispr_read_counts.command_line import cli
from crispr_read_counts.version import version

SUB_COMMANDS = ['count-single', 'merge-single', 'count-dual']

def run_command(args: List[str]):
  runner = CliRunner()
  return runner.invoke(cli, args)

@pytest.mark.parametrize('args, expected_output', [
  ([], 'Usage: cli [OPTIONS] COMMAND [ARGS]...'),
  (['--version'], f'cli, version {version}'),
  ([SUB_COMMANDS[0], '--help'], f'Usage: cli {SUB_COMMANDS[0]} [OPTIONS]'),
  ([SUB_COMMANDS[1], '--help'], f'Usage: cli {SUB_COMMANDS[1]} [OPTIONS]'),
  ([SUB_COMMANDS[2], '--help'], f'Usage: cli {SUB_COMMANDS[2]} [OPTIONS]'),
])
def test_basics(args, expected_output):
  result = run_command(args)
  assert result.exit_code == 0
  assert result.output.split('\n')[0] == expected_output
