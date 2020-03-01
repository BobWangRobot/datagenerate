#!/bin/bash
#Program:
#   generate helix
PATH=$PATH:...
export PATH
for i in {1..50};do
     phenix.python test8_4.py all 8.1  per$i.pdb 12perfect_helix.pdb
     done
