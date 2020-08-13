#!/usr/bin/env bash
set -e
pytest --cov-branch --cov-report term --cov-report html --cov=crispr_read_counts --cov-fail-under=80 -x
set +e

# these should not die:

echo -e "\n#################################"
echo      "# Running pycodestyle (style)   #"
echo      "#################################"
pycodestyle crispr_read_counts

echo -e "\n#########################################"
echo      "# Running radon (cyclomatic complexity) #"
echo      "#########################################"
radon cc -nc crispr_read_counts

echo -e "\n#########################################"
echo      "# Running radon (maintainability index) #"
echo      "#########################################"
radon mi -s crispr_read_counts

exit 0 # don't die based on assements of code quality
