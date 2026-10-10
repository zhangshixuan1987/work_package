#!/bin/sh

files="./*.png"

for file in ${files};do

filnam=`basename ${file}`

convert ${file} -trim ${filnam}_crop.png

mv ${filnam}_crop.png ${filnam}
done

#rm -rvf *crop*.pdf

#for file in fig*.pdf;do

#pdfcrop --margins '5 5 5 5' $file

#done

