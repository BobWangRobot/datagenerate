#!/bin/bash
#Program:
#   generate helix
PATH=$PATH:...
export PATH
#cd ~/Desktop/datagenerate/test/pertorb
for i in {1..50};do
     phenix.python test8_1.py back 8.1 12perfect_helix.pdb per$i.pdb
     done
