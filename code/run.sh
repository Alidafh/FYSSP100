#!/bin/bash

if [ ! -d output/figures ]; then
  mkdir -p output/figures;
  echo "created output/figures"
fi

if [ ! -d output/histograms ]; then
  mkdir -p output/histograms;
  echo "created output/histograms"
fi

FILE1=output/histograms/test_statistic_distribution_sfb1_sfs1_toys10000_bin1.root
FILE2=output/histograms/test_statistic_distribution_sfb1.5_sfs1.5_toys10000_bin10.root
FILE3=output/histograms/test_statistic_distribution_sfb2.0_sfs2.0_toys10000_bin10.root
FILE4=output/histograms/test_statistic_distribution_sfb5.0_sfs5.0_toys10000_bin10.root

if [ -f "$FILE1" -a -f "$FILE2" -a -f "$FILE3" -a -f "$FILE4" ]; then
  echo "$FILE1 exists"
  echo "$FILE2 exists"
  echo "$FILE3 exists"
  echo "$FILE4 exists"
  echo "   "
  echo "To use these distributions in the analysis type 0 (recommended)"
  echo "To calculate new distributions press any key:"
else
  echo " "
  echo "Full analysis started (perfect time for coffee...)"
  echo "--------------------------------------------------"
  python main.py 1
fi

read type

if [ $type -eq 0 ]; then
   echo "Quick run started using the existing distributions"
   echo "--------------------------------------------------"
   python main.py 0
else
  echo " "
  echo "Full analysis started (perfect time for coffee...)"
  echo "--------------------------------------------------"
  python main.py 1
fi
