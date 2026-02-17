#!/bin/sh

for file in *.png;do

filnam=`basename ${file}`

convert ${file} -trim -density 900 -bordercolor white -border 20x20 -density 900 ${filnam}_crop.png

mv ${filnam}_crop.png ${filnam}

done

