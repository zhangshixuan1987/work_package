#!/bin/sh
format="pdf"

if [ $format = "pdf" ];then

files="./*.pdf"

for file in ${files};do

filnam=`basename ${file}`

subnam=${filnam%.*}

pdfcrop $file

echo ${subnam}

mv ${subnam}-crop.pdf ${filnam}

done

fi

if [ $format = "png" ];then

files="./*.png"

for file in ${files};do

filnam=`basename ${file}`

convert ${file} -trim ${filnam}_crop.png

mv ${filnam}_crop.png ${filnam}

done

fi 


