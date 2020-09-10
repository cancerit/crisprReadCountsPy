# crisprReadCountsPy

Gets CRISPR read counts from cram files and merges readcounts from multiple read count files.

| Master | Dev |
|---|---|
|  [![Build Status](https://travis-ci.org/cancerit/crisprReadCountsPy.svg?branch=master)](https://travis-ci.org/cancerit/crisprReadCountsPy) | [![Build Status](https://travis-ci.org/cancerit/crisprReadCountsPy.svg?branch=develop)](https://travis-ci.org/cancerit/crisprReadCountsPy) |

[![Docker Repository on Quay](https://quay.io/repository/wtsicgp/crisprreadcountspy/status "Docker Repository on Quay")](https://quay.io/repository/wtsicgp/crisprreadcountspy)

## Usage

To see the full usage:

```
crisprReadCounts --help
```

## Installation

```
VERSION=X.X.X
pip install https://github.com/cancerit/crisprReadCountsPy/archive/${VERSION}.tar.gz
```

## Development environment

### Setup VirtualEnv

```
cd $PROJECTROOT
hash virtualenv || pip3 install virtualenv
virtualenv -p python3 env
source env/bin/activate
pip install -r requirements.txt
python setup.py develop # so bin scripts can find module
```

### For testing/coverage (./run_tests.sh)

```
source env/bin/activate # if not already in env
pip install -r test-requirements.txt
```

---

LICENCE
=======
Copyright (c) 2014-2020 Genome Research Ltd.

Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>

This file is part of crisprReadCounts.

crisprReadCounts is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’."
