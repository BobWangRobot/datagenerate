#! /bin/bash
cd ~/datagenerate/test13/high
function read_dir()
{
for file in `ls $2`
do
  phenix.python run_new.py $2$file
 echo $2$file
done
} 
#读取第一个参数
read_dir $1