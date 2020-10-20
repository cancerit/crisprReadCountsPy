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
  result = True
  if os.path.isfile(fn):
    result = os.access(fn, os.R_OK)
  else:
    result = False

  if not result:
    if msg_if_fail:
      sys.exit(error_msg(msg_if_fail))
    else:
      sys.exit()


def check_file_writable(fn, msg_if_fail=None):
  result = True
  if os.path.exists(fn):
    # path exists
    if os.path.isfile(fn):  # is it a file or a dir?
      # also works when file is a link and the target is writable
      result = os.access(fn, os.W_OK)
    else:
      result = False  # path is a dir, so cannot write as a file
  else:  # target does not exist, check perms on parent dir
    pdir = os.path.dirname(fn)
    if not pdir:
      pdir = '.'
    # target is creatable if parent dir is writable
    result = os.access(pdir, os.W_OK)

  if not result:
    if msg_if_fail:
      sys.exit(error_msg(msg_if_fail))
    else:
      sys.exit()
