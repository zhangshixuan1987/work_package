#!/bin/sh

rm -rvf *crop*.pdf

for file in *.pdf;do
 filnam="${file%.*}"
 if [ -f $file ];then 
   pdfcrop $file
   #pdfcrop --margins '5 5 5 5' $file
   mv $filnam-crop.pdf $filnam.pdf
 fi 
done

for file in *.png;do
 filnam="${file%.*}"
 if [ -f $file ];then
  convert -units PixelsPerInch $file -density 600 ${filnam}.png
 fi
done

