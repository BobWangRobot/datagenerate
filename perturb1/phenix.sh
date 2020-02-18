#!/bin/bash
#Program:
#   generate helix
PATH=$PATH:...
export PATH
for i in {1..50};do
     phenix.python test8_1.py all 8.1 12perfect_helix.pdb per$i.pdb
     done
