#!/bin/sh

rm -rvf *crop*.pdf

for file in *.pdf;do

filnam="${file%.*}"

pdfcrop --margins '5 5 5 5' $file

mv $filnam-crop.pdf $filnam.pdf

convert  -density 300 $filnam.pdf $filnam.png

rm -rvf $filnam.pdf

done
