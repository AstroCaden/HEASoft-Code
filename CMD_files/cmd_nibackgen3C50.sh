#!/usr/bin/env bash
set -e

nibackgen3C50 \
  "rootdir=$1" \
  "obsid=$2" \
  "bkgidxdir=${CALDB}/data/nicer/xti/bcf/bkg" \
  "bkglibdir=${CALDB}/data/nicer/xti/bcf/bkg" \
  gainepoch=2020 \
  clobber=YES
