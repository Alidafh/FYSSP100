#!/bin/bash

if [ ! -d output/figures ]; then
  mkdir -p output/figures;
  echo "created output/figures"
fi

if [ ! -d output/histograms ]; then
  mkdir -p output/histograms;
  echo "created output/histograms"
fi

FILE=output/histograms/test_statistic_distribution_sfb1_sfs1_toys10000_bin1.root

if [ -f "$FILE" ]; then
  echo "$FILE exists"
  echo "To use this type 0:"

else
  python main.py 1
fi

read type

if [ $type -eq 0 ]; then
   echo "Quick run started"
   python main.py 0
else
  python main.py 1
fi
