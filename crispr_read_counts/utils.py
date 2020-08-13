import gzip
import re
from contextlib import contextmanager
import os

PLASMID_COUNT_HEADER = re.compile(r'^sgRNA\tgene', flags=re.I)
DNA_PATTERN = re.compile(r'^[ATGC]+$', flags=re.I)
SAFE_SEQ_FORMAT = re.compile(r'^F\d+$')


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


def rev_compl(dna: str):
  trans_dict = {65: 84, 67: 71, 71: 67, 84: 65, 97: 116, 99: 103, 103: 99, 116: 97}
  return ''.join(reversed(dna)).translate(trans_dict)


def check_file_readable(fn):
  if os.path.isfile(fn):
    return os.access(fn, os.R_OK)
  else:
    return False


def check_file_writable(fn):
  '''the function is clearly borrowed from web'''
  if os.path.exists(fn):
    # path exists
    if os.path.isfile(fn):  # is it a file or a dir?
      # also works when file is a link and the target is writable
      return os.access(fn, os.W_OK)
    else:
      return False  # path is a dir, so cannot write as a file
  # target does not exist, check perms on parent dir
  pdir = os.path.dirname(fn)
  if not pdir:
    pdir = '.'
  # target is creatable if parent dir is writable
  return os.access(pdir, os.W_OK)
