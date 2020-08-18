#!/bin/bash
#Program:
#   generate helix
PATH=$PATH:...
export PATH
cd ~/Desktop/datagenerate/test4/perturb
for i in {1..50};do
     phenix.dynamics 12perfect_helix.pdb temperature=1000 number_of_step=2000 stop_at_diff=$i/10
     mv 12perfect_helix_shaken.pdb $i.pdb
     done
