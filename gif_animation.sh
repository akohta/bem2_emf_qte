#!/bin/bash

# settings
Delay=10
Loop=0

if [ $# -ne 1 ]; then
 echo "Usage : $0 image_directory"
 exit
fi

if [ ! -d $1 ]; then
 echo "$1 does not exist or is not a directory."
 exit
fi

convert -delay $Delay -loop $Loop $1/xy_Ez_*.png xy_Ez.gif &
convert -delay $Delay -loop $Loop $1/xy_Hx_*.png xy_Hx.gif &
convert -delay $Delay -loop $Loop $1/xy_Hy_*.png xy_Hy.gif &

wait
