import gzip
import re
from contextlib import contextmanager

PLASMID_COUNT_HEADER = re.compile(r'^sgRNA\tgene', flags=re.I)


def error_msg(msg: str):
  return f'#------\n# Error: {msg}\n#------'


def warning_msg(msg: str):
  return f'#------\n# Warning: {msg}\n#------'


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
