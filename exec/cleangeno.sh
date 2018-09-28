#!/usr/bin/env bash

args=("$@")

awk 'BEGIN{d=0}
{if($1!=d){
 for (x in h){
print d" "x" "h[x]};
 d=$1;
delete h;} 
h[$2]=$3" "$4; }
END{for (x in h){print $1" "x" "h[x]}}' ${args[0]}
