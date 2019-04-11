#!/usr/bin/bash

data_path=${1}
threadnum=10
tmp="/tmp/$$.fifo"
mkfifo ${tmp}
exec 6<> ${tmp}
rm ${tmp}
for((i=0; i<${threadnum}; i++))
do
    echo ""
done >&6
for file in $(ls ${data_path}/)
do
  read -u6
  {
    echo ${file}
    name="${file}/${file}.txt"
    python encode.py -f `pwd`/${data_path}/${name} -shape `pwd`/mappers/shape.5mer.txt
    
    echo "" >&6
  }&
done
wait
exec 6>&-

