#!/bin/bash
#Program:
#   generate helix
PATH=$PATH:...
export PATH
cd ~/Desktop/datagenerate/test6/per/
cp per_model.py polyGly-helix.pdb ~/Desktop/datagenerate/test6/per/structure/
cd structure/ 
for i in {1..50};do
     phenix.python per_model.py  polyGly-helix.pdb i 
     rm -rf per_model.py polyGly-helix.pdb
     cd ~/Desktop/datagenerate/test6/per/
     phenix.python index.py
     done
