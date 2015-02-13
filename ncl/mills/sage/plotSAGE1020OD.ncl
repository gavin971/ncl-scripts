load "./readSAGE1020OD.ncl"

od=readSAGE1020OD
lat=od&lat
date=od&time

a = systemfunc("mkdir -p $IMAGE_PATH/sage")
pdffile = "$IMAGE_PATH/sage/SAGE1020OD"
wks = gsn_open_wks("pdf", pdffile)

tmin=yyyymmdd_to_yyyyfrac(19910601,0)

res                     = True          ; plot mods desired

plot = gsn_csm_lat_time(wks,od(lat|:, time|:),res)
