;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;************************************************
begin


radius=(/0.4, 0.52, 0.69, 0.88, 1.2, 1.7, 2.2, 2.9, 3.7, 4.7, 5.9, \
         7.4, 9.4, 12, 15, 19, 24, 30, 37, 47, 59, 75, 94, 120, 150, \
         190, 240, 300, 380, 480, 600, 750, 950, 1200, 1500, 1900/)
radius@units="nm"
radius@long_name="Aerosol radius (nm)"

wtpct=(/28.8, 34.8, 39.8, 45.3, 47.9, 59.5, 64.7, 68.8, 70.5, 70.5, \
      70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, \
      70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, \
      70.5, 70.5, 70.5, 70.5, 70.5, 70.5/)
wtpct@long_name="Weight % H~B~2~N~O"
      
;************************************************
; plotting parameters
;************************************************
 wks   = gsn_open_wks ("pdf","kelvinH2Oplot")                ; open workstation

 res                  = True                     ; plot mods desired
; res@tiMainString     = "Weight % H2O"          ; add title
 res@trXLog = True 
 
 res@trYMinF = 25
 res@trYMaxF = 75
 
 res@trXMinF = 0.4  ; 0.1 to 100 nm range x-axis
 res@trXMaxF = 40
 
 res@tmXBMode = "Explicit"
 res@tmXBValues = (/0.4,0.6,0.8,1.,2.,4.,6.,8.,10.,20.,40./)
 res@tmXBLabels = (/"0.4",".6"," ","1","2","4","6"," ","10","20","40"/)
  res@tiXAxisString = "Aerosol Radius (nm)"
 
 res@tmEqualizeXYSizes = True
; res@tmXTLabelFont = 25 ; change the font to Times Roman
 res@tmXTLabelFontHeightF = 0.03 ; Change the font size (default 0.02)
 
; res@tiMainFont = 25 ; change the font to Times Roman
 res@tiMainFontHeightF = 0.04 ; change the font size
 res@tiUseMainAttributes = True ; use main tile font attributes for X and Y titles

 res@xyLineThicknessF = 6.0

 plot  = gsn_csm_xy (wks,radius,wtpct,res) ; create plot
end
