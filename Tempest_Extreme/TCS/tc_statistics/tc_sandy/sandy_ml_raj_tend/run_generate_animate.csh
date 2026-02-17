#!/bin/csh

cd sub_fig
sh run_crop.sh
cd ../

rm -rvf animate_sandy_utend.gif
convert -delay 30 -loop 0 sub_fig/utend_animation_tc*png  animate_sandy_utend.gif

rm -rvf animate_sandy_vtend.gif
convert -delay 30 -loop 0 sub_fig/vtend_animation_tc*png  animate_sandy_vtend.gif

