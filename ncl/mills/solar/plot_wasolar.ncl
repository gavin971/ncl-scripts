load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

datapath="$MODEL_PATH/ccsm/inputdata/atm/waccm/phot/"
imagepath="$IMAGE_PATH/solar/"

ncf = addfile(datapath+"wasolar_1948-2140_daily_c100804.nc","r")
ncf2 = addfile(datapath+"wasolar_c100106.nc","r")

date=ncf->date
date2=ncf2->date
year= date/10000.d0
print(year(0))

wks = gsn_open_wks("pdf", imagepath+"f107_1947-2140")

data=new((/2,dimsizes(date)/),double)
data(0,:)=ncf->f107
data(1,0:dimsizes(date2)-1)=ncf2->f107

  res                   = True                      ; plot mods desired
  res@trXMinF = 1940.d0
  res@trXMaxF = 2141.d0


  res@tmXBLabelFontHeightF  = 0.015
  res@gsnMaximize          = True
  res@gsnPaperOrientation = "landscape"
  res@tmXBMinorPerMajor = 10

  res@tiYAxisFontThicknessF = 0.1
;  res@tiYAxisString="TSI (W/m~S~2~N~)"
  res@xyLineColors      = (/"black","blue","red","green"/)
  res@xyDashPatterns      = (/0,0,1/)

plot1 = gsn_xy(wks, year, data, res)

;========================================

delete(data)
data=new((/2,dimsizes(date)/),double)

data(0,:)=ncf->f107a
data(1,0:dimsizes(date2)-1)=ncf2->f107a
plot2 = gsn_xy(wks, year, data, res)

;========================================

delete(data)
data=new((/2,dimsizes(date)/),double)

data(0,:)=ncf->kp
data(1,0:dimsizes(date2)-1)=ncf2->kp
plot3 = gsn_xy(wks, year, data, res)

;========================================

delete(data)
data=new((/2,dimsizes(date)/),double)

data(0,:)=ncf->ap
data(1,0:dimsizes(date2)-1)=ncf2->ap
plot4 = gsn_xy(wks, year, data, res)

;========================================

delete(data)
data=new((/2,dimsizes(date)/),double)

data(0,:)=ncf->isn
data(1,0:dimsizes(date2)-1)=ncf2->isn
plot5 = gsn_xy(wks, year, data, res)
