#!/bin/bash
#Program:
#   generate perturb
set -v
PATH=$PATH:...
export PATH
cd ~/Desktop/datagenerate/test/pertorb
for i in {1..50};do
     phenix.pdbtools 12perfect_helix.pdb shake=$i/10
     mv 12perfect_helix.pdb_modified.pdb per$i.pdb
     done
