#!/bin/bash


echo '        writed by Jiale Chen       '

echo '        start analyzing the motifs       '
echo '        5       '
echo '        4       '
echo '        3       '
echo '        2       '
echo '        1       '

mkdir -p 
ls *_[12].bed |sed 's/_[12].bed//' |sort -u |while read line ;do 
	echo  '        making new dictory'
	mkdir -p ../motif/${line}_1/
	
	findMotifsGenome.pl ${line}_1.bed hg38 ../motif/${line}_1/ -size 200 -mask 
	
	echo  '        making new dictory'
	mkdir -p ../motif/${line}_2/
	
 	findMotifsGenome.pl ${line}_2.bed hg38 ../motif/${line}_2/ -size 200 -mask 
  
done 
echo  '   finished   '
echo  ` date `
