# Changes

## 2.0.1

* Used the Quay.io repo for original Perl code for this code base 

## 2.0.0

* migrated Perl code of crisprReadCounts to Python3
* Comparing to Perl version: sub-command "single-count" will not count non-primary alignment and vendor failed alignments
* Comparing to Perl version: sub-command "merge-single" option `--plasmid, -p` is a flag, requiring NO parameter.
* Added an option to set library file delimiter for `count-single` to use. By default, it assumes library file is tab delimited.
* Sub command can optionally generate a stats file in JSON for `count-single`.
* Added unit tests and some code improvements.
