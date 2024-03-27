#!/bin/bash

mkdir -p 
ls *_[12].bed |sed 's/_[12].bed//' |sort -u |while read line ;do 
	mkdir -p ../motif/${line}_1/
	findMotifsGenome.pl ${line}_1.bed hg38 ../motif/${line}_1/ -size 200 -mask 
	mkdir -p ../motif/${line}_2/
 	findMotifsGenome.pl ${line}_2.bed hg38 ../motif/${line}_2/ -size 200 -mask 
  
done 
