load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "../calendar_decode2.ncl"

path="$MODEL_PATH/ccsm/inputdata/atm/cam/solar/"
a=addfile(path+"spectral_irradiance_Lean_1950-2140_daily_c100804.nc","r")

time=a->time
date=calendar_decode2(time,4)
data=a->tsi

filename="$IMAGE_PATH/WA4solarPlot"

  wks   = gsn_open_wks ("pdf",filename) ; open workstation

  res                   = True                      ; plot mods desired
  res@xyLineThicknessF = 1.0
  res@gsnPaperOrientation = "landscape"
  res@gsnMaximize          = True
;  res@vpHeightF             = 0.5
;  res@vpWidthF              = 1.2

  res@xyLineColor       = ("red")          ; change line color
  res@xyMonoLineColor   = True            ; same color for all lines
  res@xyDashPattern     = 0                         ; all solid lines
  res@trYMinF = 1357.5
  res@trYMaxF = 1363.5

  res@trXMinF = 2010
  res@trXMaxF = 2090

  res@tiYAxisString = "Total Solar Irradiance (W/m~S~2~N~)"
  res@tiXAxisString = "Year"
;  res@tmXTLabelFontHeightF = 0.04 ; Change the font size (default 0.02)

;    res@tiMainFont = 25 ; change the font to Times Roman
;    res@tiMainFontHeightF = 0.04 ; change the font size
;    res@tiUseMainAttributes = True ; use main tile font attributes for X and Y titles


  plot = gsn_csm_xy (wks,date,data,res) ; draw second plot
  
  delete(plot)
  delete(wks)
  
print("mutt -a " + filename + ".pdf mikejmills@mac.com < /dev/null")
mutter=systemfunc("mutt -a " + filename + ".pdf mikejmills@mac.com < /dev/null")
print(mutter)
