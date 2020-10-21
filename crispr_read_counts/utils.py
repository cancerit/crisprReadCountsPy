import gzip
import re
from contextlib import contextmanager
import os
import sys

PLASMID_COUNT_HEADER = re.compile(r'^sgRNA\tgene', flags=re.I)
DNA_PATTERN = re.compile(r'^[ATGC]+$', flags=re.I)
SAFE_SEQ_FORMAT = re.compile(r'^F\d+$')

dna_complement_tr_table = str.maketrans('ACGTacgt', 'TGCAtgca')


def rev_compl(dna: str) -> str:
    return dna[::-1].translate(dna_complement_tr_table)


def error_msg(msg: str):
  return f'#------\n# Error: {process_multiple_lines(msg)}\n#------'


def warning_msg(msg: str):
  return f'#------\n# Warning: {process_multiple_lines(msg)}\n#------'


def process_multiple_lines(msg):
  return "\n".join(
    [
      line if index == 0 else '# ' + line
      for index, line in enumerate(msg.split('\n'))
    ]
  )


@contextmanager
def open_plain_or_gzipped_file(file: str):
    if file.endswith('.gz'):
      f = gzip.open(file, 'rt')
    else:
      f = open(file, 'r')
    try:
        yield f
    finally:
        f.close()


def check_file_readable(fn, msg_if_fail=None):
  result: bool = os.path.isfile(fn) and os.access(fn, os.R_OK)

  if not result:
    if msg_if_fail:
      sys.exit(error_msg(msg_if_fail))
    else:
      sys.exit()


def check_file_writable(fn, msg_if_fail=None):
  result: bool = (
    (os.path.isfile(fn) and os.access(fn, os.W_OK)) if os.path.exists(fn) else
    os.access(os.path.dirname(fn) or '.', os.W_OK)
  )

  if not result:
    if msg_if_fail:
      sys.exit(error_msg(msg_if_fail))
    else:
      sys.exit()
