#!/bin/csh

set stim = "2016012000"
set etim = "2016012418"

cd trk_sub_fig
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

rm -rvf ../animate_etc_trk*.gif
echo $pfile 

convert -delay 50 -loop 0 $pfile ../animate_etc_trk_2016.gif
