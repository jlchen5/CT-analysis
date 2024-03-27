#!/bin/bash

mkdir -p 
ls *_[12].bed |sed 's/_[12].bed//' |sort -u |while read line ;do 
	mkdir -p ../motif/$line/
	findMotifsGenome.pl ${line}_1.bed hg38 ../motif/$line/ -size 200 -mask 
done 
