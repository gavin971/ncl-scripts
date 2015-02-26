load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
  
ostype=getenv("OSTYPE")
host=getenv("HOSTNAME")

sagepath="/glade/p/work/mmills/data/SAGE_II_v7.00/CCMI_Filled_Data_V2/netCDF/"
if (ostype.eq."darwin") then
  sagepath="/Volumes/Data/Data/SAGE2/SAGE_II_v7.00/CCMI_Filled_Data_V2/netCDF/"
end if

a=addfile(sagepath+"ODstrat.nc","r")
ODst1020=a->ODst1020
ODst525=a->ODst525
ODst452=a->ODst452
ODst386=a->ODst386

lat=ODst1020&lat
date=ODst1020&time

print("mkdir -p $IMAGE_PATH/sage")
b = systemfunc("mkdir -p $IMAGE_PATH/sage")
pdffile = "$IMAGE_PATH/sage/SAGE_OD"
wks = gsn_open_wks("pdf", pdffile)

res                     = True          ; plot mods desired
;res@trXMinF=tmin
res@trXMinF=1991.0
res@trXMaxF=2000.0
res@gsnPaperOrientation = "landscape"
res@vpWidthF            = 1.2           ; change aspect ratio of plot
res@vpHeightF           = 0.55
res@cnFillOn            = True          ; turn on color
res@gsnSpreadColors     = True          ; use full range of colormap
res@gsnMaximize         = True          ; expand plot size to maximum
res@cnLevelSelectionMode = "ExplicitLevels" 
res@cnLevels = 10.^fspan(-.5,-3.,11)

;res@gsnCenterString=ODst1020@long_name
plot = gsn_csm_lat_time(wks,ODst1020(lat|:, time|:),res)
delete(plot)

;res@gsnCenterString=ODst525@long_name
plot = gsn_csm_lat_time(wks,ODst525(lat|:, time|:),res)
delete(plot)

;res@gsnCenterString=ODst452@long_name
plot = gsn_csm_lat_time(wks,ODst452(lat|:, time|:),res)
delete(plot)

;res@gsnCenterString=ODst386@long_name
plot = gsn_csm_lat_time(wks,ODst386(lat|:, time|:),res)
delete(plot)
