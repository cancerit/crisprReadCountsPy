import click
from .version import version
import sys


@click.group()
@click.version_option(version)
def cli():
  pass


@cli.command()
@click.option(
  # TODO: it should not be "dir", but perl code used the word, change it in later version.
  '--dir', '-i',
  metavar='FILE',
  required=True,
  help='Input sample CRAM file.')
@click.option(
  '--plasmid', '-p',
  metavar='FILE',
  help='Plasmid count tsv file.')
@click.option(
  '--library', '-l',
  metavar='FILE',
  required=True,
  help='Input sample CRAM file.')
@click.option(
  '--output', '-o',
  metavar='FILE',
  required=True,
  help='Output read counts file.')
@click.option(
  '--ref', '-r',
  metavar='FILE',
  required=True,
  help='Genome reference FASTA (e.g.: genome.fa) file.')
@click.option(
  '--trim', '-t',
  metavar='INT',
  default=0,
  help='Remove N bases of leading sequence.')
@click.option(
  '--reverse-complement', '-rc',
  default=False,
  is_flag=True,
  help='Reverse complementing reads when mapping to guide sequences (reads are reverse complemented prior to trimming).')
def count_single(**kwargs):
  from .count import count_single
  count_single(kwargs)


@cli.command()
@click.option(
  '--library', '-l',
  metavar='FILE',
  required=True,
  help='Guide library file.')
@click.option(
  '--fastq1', '-f1',
  required=True,
  metavar='FILE',
  help='R1 fastq file.')
@click.option(
  '--fastq2', '-f2',
  required=True,
  metavar='FILE',
  help='R2 fastq file.')
@click.option(
  '--sample', '-n',
  metavar='STRING',
  required=True,
  help='Sample name.')
@click.option(
  '--reads', '-r',
  metavar='FILE',
  required=True,
  help='Output classified reads file.')
@click.option(
  '--stats', '-s',
  metavar='FILE',
  required=True,
  help='Output classification stats file.')
@click.option(
  '--counts', '-c',
  metavar='FILE',
  required=True,
  help='Output read counts result file.')
def count_dual(**kwargs):
  from .count import count_dual
  count_dual(kwargs)


@cli.command()
@click.option(
  '--input', '-i',
  metavar='FILE',
  # NOTE: not the native way to handle list of input. It's to repect the existing Perl version interface.
  help='Comma separated list of input files.')
@click.option(
  '--output', '-o',
  metavar='FILE',
  help='Output file.')
@click.option(
  '--plasmid', '-p',
  metavar='FILE',
  default=False,
  is_flag=True,
  help='Has plasmid counts.')
def merge_single(**kwargs):
  from .merge import merge_single
  merge_single(kwargs)


def main():
  cli()
