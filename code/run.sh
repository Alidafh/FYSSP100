#!/bin/bash

if [ ! -d output/figures ]; then
  mkdir -p output/figures;
  echo "created output/figures"
fi

if [ ! -d output/histograms ]; then
  mkdir -p output/histograms;
  echo "created output/histograms"
fi

echo "for quick run type 0"

read type

if [ $type -eq 0 ]; then
   echo "quick run started"
   python main.py 0

fi
