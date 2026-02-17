#!/bin/csh

cd sub_fig
sh run_crop.sh
cd ../

rm -rvf animate*.gif
convert -delay 30 -loop 0 sub_fig/animation_tc*png  animate_sandy.gif
