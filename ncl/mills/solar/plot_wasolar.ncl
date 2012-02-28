load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

datapath="$MODEL_PATH/ccsm/inputdata/atm/waccm/solar/"
imagepath="$IMAGE_PATH/solar/"

ncf2 = addfile(datapath+"wasolar_1948-2140_daily_c100805.nc","r")
ncf  = addfile(datapath+"wasolar_1947-2011_daily_c120227.nc","r")

date=ncf->date
date2=ncf2->date

dayFirst=20110101
dayLast=date(dimsizes(date)-1)
print("dayLast="+dayLast)
dayInd=ind(date.ge.dayFirst.and.date.le.dayLast)
ndays=dimsizes(dayInd)
print("ndays="+ndays)
dayInd2=ind(date2.ge.dayFirst.and.date2.le.dayLast)
ndays2=dimsizes(dayInd2)
print("ndays2="+ndays2)

year = yyyymmdd_to_yyyyfrac(date(dayInd), 0.)
print("year="+year(0)+"-"+year(ndays-1))

wks = gsn_open_wks("pdf", imagepath+"f107_2011")

data=new((/2,ndays/),double)
data(0,:)=ncf->f107(dayInd)
data(1,:)=ncf2->f107(dayInd2)

  res                   = True                      ; plot mods desired
  res@trXMinF = 2011.d0
  res@trXMaxF = 2012.d0


  res@tmXBLabelFontHeightF  = 0.015
  res@gsnMaximize          = True
  res@gsnPaperOrientation = "landscape"
  res@tmXBMinorPerMajor = 10

  res@tiYAxisFontThicknessF = 0.1
;  res@tiYAxisString="TSI (W/m~S~2~N~)"
  res@xyLineColors      = (/"black","red"/)
  res@xyDashPatterns      = (/0,1/)

plot1 = gsn_xy(wks, year, data, res)

;========================================

delete(data)
data=new((/2,ndays/),double)

data(0,:)=ncf->f107a(dayInd)
data(1,:)=ncf2->f107a(dayInd2)
plot2 = gsn_xy(wks, year, data, res)

;========================================

delete(data)
data=new((/2,ndays/),double)

data(0,:)=ncf->kp(dayInd)
data(1,:)=ncf2->kp(dayInd2)
plot3 = gsn_xy(wks, year, data, res)

;========================================

delete(data)
data=new((/2,ndays/),double)

data(0,:)=ncf->ap(dayInd)
data(1,:)=ncf2->ap(dayInd2)
plot4 = gsn_xy(wks, year, data, res)

;========================================

delete(data)
data=new((/2,ndays/),double)

data(0,:)=ncf->isn(dayInd)
data(1,:)=ncf2->isn(dayInd2)
plot5 = gsn_xy(wks, year, data, res)
