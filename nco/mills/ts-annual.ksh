#!/bin/ksh

;runname='b40.1850.track1.2deg.wcm.007'
runname='b40.1850.2deg.wset.003'
wadir='/data3/mmills/'$runname'/h0'
anndir='/data3/mmills/'$runname'/h0.annual'
searchstr='cam2.h0'
addstr='TS.annual'

year1=1
year2=18

# script to average and standard deviation of waccm h0 files

mkdir -p $anndir

year=$year1
while [[ $year -le $year2 ]];do

  if [[ year -lt 10 ]]
    then
      yearstr='000'$year
    elif [[ year -lt 100 ]]
    then
      yearstr='00'$year
    elif [[ year -lt 1000 ]]
    then
      yearstr='0'$year
    else
      yearstr=$year
    fi
  
  ncfiles=$wadir'/'$runname'.'$searchstr'.'$yearstr*'.nc'
  avgfile=$anndir'/'$runname'.'$searchstr'.'$yearstr'.'$addstr'.nc'

  # create annual mean
  echo averaging $ncfiles
  ncea -v TS,date $ncfiles $avgfile
  echo created $avgfile
   
  (( year+=1 ))
     
done

echo "Annual averaging complete.  Now run fixdate.ncl to fix the dates!"

export runname
export anndir
export year1
export year2
export searchstr
export addstr

ncl ~/ncl-scripts/fixdate.ncl

#catfile=$anndir'/'$runname'.'$searchstr'.'$year1'-'$year2'.TS.nc'
#echo creating $catfile
#ncrcat -O $anndir'/*.'$addstr'.nc' $catfile
