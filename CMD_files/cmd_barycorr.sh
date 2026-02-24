#!/usr/bin/env bash
set -e


export HEADASNOQUERY=1
export HEADASYES=1

barycorr \
  "infile=$1" \
  "outfile=$2" \
  "orbitfiles=$3" \
  "ra=$4" \
  "dec=$5" \
  refframe=ICRS \
  ephem=JPLEPH.430 \
  "barytime=$6" \
  clobber=yes
