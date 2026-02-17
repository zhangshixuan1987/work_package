#!/bin/sh

files="./*.png"

for file in ${files};do

filnam=`basename ${file}`

convert ${file} -trim ${filnam}_crop.png

mv ${filnam}_crop.png ${filnam}
done

