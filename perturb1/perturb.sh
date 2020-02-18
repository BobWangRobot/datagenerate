#!/bin/bash
#Program:
#   generate helix
PATH=$PATH:...
export PATH
for i in {1..50};do
     phenix.dynamics 12perfect_helix.pdb temperature=1000 number_of_steps=500 stop_at_diff=0.6 random_seed=$i*10240
     mv 12perfect_helix_shaken.pdb per$i.pdb
     done
