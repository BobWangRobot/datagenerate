#!/bin/bash
#Program:
#   generate helix
PATH=$PATH:...
export PATH
cd ~/Desktop/datagenerate/test
read -p "please input filename:" filename
read -p "please input resseq:" seq
phenix.pdbtools $filename keep="chain A and resseq $seq"
