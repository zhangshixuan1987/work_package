#!/bin/csh

set stim = "2017020300"
set etim = "2017021112"

cd prec_sub_fig
cp ../run_crop.sh .
sh run_crop.sh

set files = `ls *.png`
set nfile = $#files

set ftim =
set pfile =
set i = 1
while ( $i <= $nfile)
 set tim  = `echo $files[$i] | awk '{print substr($0,10,10)}'`
 set ftim = ($ftim $tim)

 if ( $tim >= $stim && $tim <= $etim ) then
  set pfile = ($pfile $files[$i])
 endif

 @ i++
end

rm -rvf ../animate_prec*.gif
convert -delay 30 -loop 0 $pfile ../animate_prec_201702.gif

