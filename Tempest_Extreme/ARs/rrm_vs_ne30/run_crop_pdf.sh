#!/bin/sh

for file in *.png;do

filnam=`basename ${file}`

convert ${file} -trim -bordercolor white -border 20x20 -density 1200 ${filnam}_crop.png

mv ${filnam}_crop.png ${filnam}

done

