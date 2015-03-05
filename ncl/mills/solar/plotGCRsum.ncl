load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

rootpath="/glade/p/cesmdata/cseg/inputdata/atm/waccm/gcrs/"
filename="gcr_prod_NO_c130920.nc"

a=addfile(rootpath+filename,"r")
prod=a->prod ; (time, lev, lat, lon)
date=a->date
yearint=date/10000
year=yearint;+0.5
year2=ispan(1955,2101,1);+0.5
prodSum2=new(dimsizes(year2),"float")
prodSum2!0="year"
prodSum2&year=year2

rad = 4.0*atan(1.0)/180.0
lat = prod&lat
gw = cos(lat*rad)

prodSum=dim_sum(wgt_areaave(prod(:, :, :, :), gw(:), 1.0, 0))
prodSum!0="year"
prodSum&year=year

print(year+" "+prodSum)

prodSum2({1960:2010})=(/prodSum({1960:2010})/)
prodSum2({1955:1959})=(/prodSum({1997:2001})/)
prodSum2({2011:2054})=(/prodSum({1967:2010})/)
prodSum2({2055:2098})=(/prodSum({1967:2010})/)
prodSum2({2099:2101})=(/prodSum({1967:1969})/)

wks = gsn_open_wks("pdf", "$IMAGE_PATH/GCRsum")

res = True                      ; plot mods desired
res@gsnDraw           = False
res@gsnFrame          = False
res@vpHeightF            = 0.5              ; change aspect ratio
res@vpWidthF             = 1.2              ; change aspect ratio
res@gsnMaximize          = True
res@gsnPaperOrientation = "landscape"
res@xyLineThicknessF=5.0
;res@trXMinF = 1960.0
;res@trXMaxF = 2011.0
res@trXMinF = 1955.0
res@trXMaxF = 2101.0
;res@tmXBMinorPerMajor = 9
;res@tmXTMinorPerMajor = 9
year3=year2+0.5

plot1 = gsn_xy(wks, year3, prodSum2, res)

res2=res
res2@xyLineColor="red"
res2@xyDashPattern=2
year4=year+0.5
overlay(plot1,gsn_xy(wks, year4, prodSum, res2))

draw(plot1)

delete(wks)

print("evince $IMAGE_PATH/GCRsum.pdf")
