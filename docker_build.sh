#! /bin/bash

set -xue

# assume that the version string is withon single quotes in the perl module file
version_string=$(grep 'version=' crispr_read_counts/version.py | cut -d "'" -f 2)

if [[ -z "${version_string}" ]]; then echo "Failed to capture version string"; exit 1; fi

docker build --build-arg VERSION=$version_string $@
