#!/bin/csh

set stim = "2017020500"
set etim = "2017021118"

cd iwv_sub_fig
cp ../run_crop.sh .
sh run_crop.sh

set files = `ls *.png`
set nfile = $#files

set ftim = 
set pfile = 
set i = 1
while ( $i <= $nfile)
 set tim  = `echo $files[$i] | awk '{print substr($0,9,10)}'`
 set ftim = ($ftim $tim)

 if ( $tim >= $stim && $tim <= $etim ) then 
  set pfile = ($pfile $files[$i]) 
 endif 

 @ i++
end 

rm -rvf ../animate_iwv*.gif
convert -delay 30 -loop 0 $pfile ../animate_iwv_201702.gif

