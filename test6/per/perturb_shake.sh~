#!/bin/bash
#Program:
#   generate helix
PATH=$PATH:...
export PATH
cd ~/Desktop/datagenerate/test4/perturb
for i in {1..50};do
     phenix.pdbtools 12perfect_helix.pdb shake=$i/10
     mv 12perfect_helix.pdb_modified.pdb per$i.pdb
     done
