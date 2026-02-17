#!/bin/sh

for file in *.png;do

filnam=`basename ${file}`

convert ${file} -trim ${filnam}_crop.png

mv ${filnam}_crop.png ${filnam}
done

