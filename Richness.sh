#!/bin/bash

Rscript ~/scripts/multiplex/richness.R

for f in *_multi_final.txt
 do NAME=${f%%_*} && T=${f%%.*} && INDEX=$(echo "$T".txt | grep -oP '(?<=_)\d+(?=\_)')
 echo $NAME;
 echo $INDEX;
 echo $T;
 awk -F, -v s="$NAME" '{$2="\t"s;print}' Richness.txt > Richness-temp.txt
 tr -d '-' < Richness-temp.txt > Richness-tempII.txt
 awk -F, -v s="$INDEX" '{$3="\t"s;print}' Richness-tempII.txt > Richness-tempIII.txt
 rm Richness.txt Richness-temp.txt Richness-tempII.txt
 mv Richness-tempIII.txt Richness.txt
 done