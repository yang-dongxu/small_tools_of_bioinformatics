#!/bin/bash

process_1(){
filename=$1
f=$(basename $1)
project=${f%%_*}
region_tmp=${f##*_}
region=${region_tmp%.*}

cat $filename | sed 's/^/'$project'\t'$region'\t/'

#echo "$filename $project $region"
}



if [ $# -gt 0 ];then
    exec 0<$1;    #将文件绑定到标准输入（0-标准输入 1-标准输出 3-标准错误），默认第一个参数是输入的文件；
fi

while read line
do
    process_1 $line
done   #从标准输入读取数据
exec 0<&-  
