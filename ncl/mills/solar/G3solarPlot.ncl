load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

path="$MODEL_PATH/ccsm/inputdata/atm/cam/solar/"
a=addfile(path+"SOLAR_TSI_Lean_1610-2140_annual_c100301.nc","r")
b=addfile(path+"SOLAR_TSI_Lean_1610-2140_annual_G3solar_c110222.nc","r")

date=a->date/1e4
data=new((/2,531/),double)
data(0,:)=a->tsi
data(1,:)=b->tsi

  wks   = gsn_open_wks ("pdf","$IMAGE_PATH/G3solarPlot") ; open workstation

  res                   = True                      ; plot mods desired
  res@xyLineThicknessF = 2.0

  res@xyLineColors       = (/"red","blue"/)          ; change line color
  res@xyMonoLineColor   = False            ; same color for all lines
  res@xyDashPattern     = 0                         ; all solid lines
  res@trYMinF = 1350
  res@trYMaxF = 1364

  res@trXMinF = 2010
  res@trXMaxF = 2090

  res@tiYAxisString = "Total Solar Irradiance (W/m~S~2~N~)"
  res@tiXAxisString = "Year"
;  res@tmXTLabelFontHeightF = 0.04 ; Change the font size (default 0.02)

;    res@tiMainFont = 25 ; change the font to Times Roman
;    res@tiMainFontHeightF = 0.04 ; change the font size
;    res@tiUseMainAttributes = True ; use main tile font attributes for X and Y titles


  plot = gsn_csm_xy (wks,date,data,res) ; draw second plot
