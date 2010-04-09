sootfile=addfile("/data5/mmills/vgsoot08/h0zm/vgsoot08.cam2.h0zm.annual.nc","r")
controlfile=addfile("/data5/mmills/vgcontrol02/h0zm/vgcontrol02.cam2.h0zm.annual.nc","r")

lat=round(sootfile->lat,0)
O3_col_soot=sootfile->O3_col
O3_col_control=controlfile->O3_col
O3_col_diff=O3_col_soot-O3_col_control
O3_col_pct_diff=100.0*O3_col_diff/O3_col_control
do i = dimsizes(lat)-1, 0, 1
  print(sprintf("%3.0f",lat(i))    +"  "+sprintf("%5.2f",O3_col_control(i))+"  "+\
    sprintf("%4.2f",O3_col_soot(i))+"  "+sprintf("%6.3f",O3_col_diff(i))   +"  "+\
    sprintf("%6.3f",O3_col_pct_diff(i)))
end do
