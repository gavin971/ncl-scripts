load "./readSAGE1020OD.ncl"

od=readSAGE1020OD
lat=od&lat
date=od&time

a = systemfunc("mkdir -p $IMAGE_PATH/sage")
pdffile = "$IMAGE_PATH/sage/SAGE1020OD"
wks = gsn_open_wks("pdf", pdffile)

tmin=yyyymmdd_to_yyyyfrac(19910601,0)

res                     = True          ; plot mods desired
;res@trXMinF=tmin
res@trXMinF=1991.0
res@trXMaxF=1995.0
res@gsnPaperOrientation = "landscape"
res@vpWidthF            = 1.2           ; change aspect ratio of plot
res@vpHeightF           = 0.55
res@cnFillOn            = True          ; turn on color
res@gsnSpreadColors     = True          ; use full range of colormap
res@gsnMaximize         = True          ; expand plot size to maximum
res@cnLevelSelectionMode = "ExplicitLevels" 
res@cnLevels = 10.^fspan(-.5,-3.,11)

plot = gsn_csm_lat_time(wks,od(lat|:, time|:),res)
