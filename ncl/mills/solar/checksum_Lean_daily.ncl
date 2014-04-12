; Program checks solar lean netCDF file that integral of spectrum
; equals TSI, and that TSI is same as for file
; spectral_irradiance_Lean_1950-2140_daily_c100803.nc

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
host   = systemfunc("hostname")

if (host.eq."modeling1") then
  datapath="/data4/mmills/inputdata/solar/"
  inputdatapath="/data4/mmills/inputdata/atm/cam/solar/"
else
  inputdatapath="/glade/p/cesmdata/cseg/inputdata/atm/cam/solar/"
  ;datapath="$MODEL_PATH/ccsm/inputdata/atm/cam/solar/"
  datapath="/glade/p/work/mmills/inputdata/solar/"
end if

imagepath="$IMAGE_PATH/solar/"
dum=systemfunc("mkdir -p "+imagepath)

;ncf = addfile(datapath+"spectral_irradiance_Lean_1978-2013_daily_GOME-Mg_Leap_c140227.nc","r")
ncf = addfile(datapath+"spectral_irradiance_Lean_1978-2014_daily_GOME-Mg_Leap_c140412.nc","r")

print("reading ncf->ssi")
ssi=ncf->ssi
print("reading ncf->tsi")
tsi=ncf->tsi
print("dimsizes(tsi)="+dimsizes(tsi))
print("min="+min(tsi))
print("max="+max(tsi))
print("reading ncf->band_width")
band_width=ncf->band_width
;print("reading ncf->time")
;time=ncf->time
print("reading ncf->date")
date=ncf->date
print("min="+min(date))
print("max="+max(date))
ndate=dimsizes(date)

print("calculating tsi_sum")
tsi_sum = tsi ; get metadata
tsi_sum = ssi # band_width
print("converting mW to W")
tsi_sum = tsi_sum/1000.d0 ; Convert mW to W
print("calculating tsi_diff")
tsi_diff=tsi-tsi_sum

;print("calculating tsi_sum2")
;tsi_sum2 = ssi(:,1:) # band_width(1:)
;print("converting mW to W")
;tsi_sum2 = tsi_sum2/1000.d0 ; Convert mW to W
;tsi_sum2!0="time"
;tsi_sum2&time=time
;print("calculating tsi_diff3")
;tsi_diff3=tsi-tsi_sum2


print("calculating year")
year = yyyymmdd_to_yyyyfrac(date, 0.)

filename=imagepath+"ssi_checksum_1978-2014"
print("plotting to "+filename+".pdf")
wks = gsn_open_wks("pdf", filename)
  res                   = True                      ; plot mods desired
  res@trXMinF = 1978.8d0
  res@trXMaxF = 2014.3d0

  res@trYMinF = 1360.2
  res@trYMaxF = 1362.2
  res@tmXBLabelFontHeightF  = 0.015
  res@gsnMaximize          = True
  res@gsnPaperOrientation = "landscape"
  res@tmXBMinorPerMajor = 10

;  res@tmYBLabelFontHeightF  = 0.01
  res@tiYAxisFontThicknessF = 0.1
;  res@tiXAxisString="Time (years)"
  res@tiYAxisString="TSI (W/m~S~2~N~)"
  res@xyLineColors      = (/"blue","red"/)
  res@xyDashPatterns      = (/0,1/)
  
  data=new((/2,ndate/),double)
  data(1,:)=tsi(:)
  data(0,:)=tsi_sum(:)
  plot1 = gsn_xy(wks, year, data, res)
  
  delete(res@trYMinF)
  delete(res@trYMaxF)
  delete(res@xyLineColors)
  delete(res@xyDashPatterns)  
  res@tiYAxisString="SSI sum/TSI"
  data2 = tsi_sum(:)/tsi(:)
  plot2 = gsn_xy(wks, year, data2, res)
  
delete(wks)
delete(res)
  
