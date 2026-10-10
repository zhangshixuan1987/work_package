#!/bin/sh

files="./*.png"

for file in ${files};do

filnam=`basename ${file}`

convert ${file} -trim ${filnam}_crop.png

mv ${filnam}_crop.png ${filnam}
done

exit

rm -rvf *crop*.pdf

for file in *.pdf;do

filnam="${file%.*}"

pdfcrop --margins '5 5 5 5' $file

mv $filnam-crop.pdf $filnam.pdf

done
