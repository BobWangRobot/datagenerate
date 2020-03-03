#!/bin/bash
#Program:
#   generate helix
PATH=$PATH:...
export PATH
for i in {1..50};do
     phenix.dynamics 12perfect_helix.pdb temperature=1000 number_of_steps=500 stop_at_diff=$i/10 random_seed=93227561
     mv 12perfect_helix_shaken.pdb per$i.pdb
     done
