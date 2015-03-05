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
year=yearint+0.5

rad = 4.0*atan(1.0)/180.0
lat = prod&lat
gw = cos(lat*rad)

prodSum=dim_sum(wgt_areaave(prod(:, :, :, :), gw(:), 1.0, 0))

print(year+" "+prodSum)

wks = gsn_open_wks("pdf", "$IMAGE_PATH/GCRsum")

res = True                      ; plot mods desired
res@vpHeightF            = 0.5              ; change aspect ratio
res@vpWidthF             = 1.2              ; change aspect ratio
res@gsnMaximize          = True
res@gsnPaperOrientation = "landscape"
res@xyLineThicknessF=5.0
res@trXMinF = 1960.0
res@trXMaxF = 2011.0
res@tmXBMinorPerMajor = 10
res@tmXTMinorPerMajor = 10

plot1 = gsn_xy(wks, year, prodSum, res)

delete(wks)

print("evince $IMAGE_PATH/GCRsum.pdf")
