; Dust-Sulfates
;
; Monthly average plots

load "tests/common.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("monave_za")
procedure monave_za()

local files, nFiles, historyFile

begin
;  runname="Solar Max"
  runname="Solar Min"
  
  print("")
  print(runname)
  print(runname .eq. "Solar Min")
  print("  Plot monthly average mixed sulfate densities ...")
  
  ; Open the  history files
  
  if (runname .eq. "Solar Min") then
  
    historyPattern="waccm/dust-sulfur/run/125smin/daily/all.125smin.cam2.h3.1998-0[67]*"
    
    historyPattern2="waccm/dust-sulfur/run/125smin/daily/all.125smin.cam2.h3.1998-01*"
    historyNames2 = systemfunc("ls " + historyPattern2)
    
  else
  
    historyPattern="waccm/dust-sulfur/run/125expuT15/daily/all.125expuT15.cam2.h3.1997-0[67]*"
    
    historyNames2="waccm/dust-sulfur/run/125expuT15b/daily/all.125expuT15b.cam2.h3.199" + \
       (/ "7-12-31-00000.nc", "8-01-05-00000.nc", "8-01-10-00000.nc", "8-01-15-00000.nc", \
          "8-01-20-00000.nc", "8-01-25-00000.nc", "8-01-30-00000.nc" /)
    
  end if
  
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  print("historyNames:"+dimsizes(historyNames))
  print(historyNames)
   
  history2 = addfiles(historyNames2, "r")
  
;  historyNames2 = historyNames
;  history2 = history
  
  ListSetType(history2, "cat")
  print("historyNames2:"+dimsizes(historyNames2))
  print(historyNames2)

 ; Compute the pressure levels ...
  time     = addfiles_GetVar(history, historyNames, "time")
;  time2    = addfiles_GetVar(history2, historyNames2, "time")
  time2   = time

  lev     = history[4]->lev
  nLevs   = dimsizes(lev)
  lat     = history[0]->lat
  nLats   = dimsizes(lat)
  lon     = history[0]->lon
  nLons   = dimsizes(lon)
  
  ; Determine an approximate dz/dp.
  hyam 	= history[0]->hyam
  hybm 	= history[0]->hybm
  hyai 	= history[0]->hyai
  hybi 	= history[0]->hybi
  P0 	  = history[0]->P0
  PS 	  = history[0]->PS
  
  T    = addfiles_GetVar(history, historyNames, "T")
  T2   = addfiles_GetVar(history2, historyNames2, "T")
  date = addfiles_GetVar(history, historyNames, "date")
  date2= addfiles_GetVar(history2, historyNames2, "date")
  
  mixsul10A  = addfiles_GetVar(history,  historyNames,  "MIXSUL10A")
  pursul10A  = addfiles_GetVar(history,  historyNames,  "PURSUL10A")
  totsul10A  = mixsul10A ; get dim wrappers
  totsul10A  = totsul10A + pursul10A
  
  mixsul10A2 = addfiles_GetVar(history2, historyNames2, "MIXSUL10A")
  pursul10A2 = addfiles_GetVar(history2, historyNames2, "PURSUL10A")
  totsul10A2 = mixsul10A2 ; get dim wrappers
  totsul10A2 = totsul10A2 + pursul10A2
  
;  time0701=ind(date.eq."19970701")
;  time0731=ind(date.eq."19970731")
  time0701=ind(date.eq."19980701")
  time0731=ind(date.eq."19980731")
 
  time0101=ind(date2.eq."19980101")
  time0131=ind(date2.eq."19980131")
  
;  time0101=time0701
;  time0131=time0731

  print("time0101:"+time0101)
  print("time0131:"+time0131)

  lev_p   = lev
  lev_p@units = "mb"
  lev_p@long_name = "Pressure (mb)"
  
  ; Make three sets of data:
  ;  control
  ;  test
  ;  difference
  month = (/ 1, 7 /)
;  month = (/ 7, 7 /)
  nMonths = dimsizes(month)
  runs = (/ "January", "July" /)
;  runs = (/ "Pure Sulfates", "Mixed Sulfates" /)
  nRuns = dimsizes(runs)
  
  Ts = new((/ nRuns, nLevs, nLats /), float)
  Ts!0 = "run"
  Ts&run = runs
  Ts!1 = "lev"
  Ts&lev = lev_p
  Ts!2 = "lat"
  Ts&lat = lat
  Ts@long_name = "Temperature"
  Ts@units = "K"
  
 ; Calculate the zonal monthly average of temperature
  print("Averaging July temperature ...")
  T_month=dim_avg_Wrap(T(lev|:,lat|:,lon|:,time|time0701:time0731))
  print("Zonally averaging July temperature ...")
  Ts(lev|:,lat|:, run|1) = dim_avg(T_month(lev|:,lat|:,lon|:))
  
  print("Averaging January temperature ...")
  T_month=dim_avg_Wrap(T2(lev|:,lat|:,lon|:,time|time0101:time0131))
  print("Zonally averaging January temperature ...")
  Ts(lev|:,lat|:, run|0) = dim_avg(T_month(lev|:,lat|:,lon|:))
  
  N1s = new((/ nRuns, nLevs, nLats /), float)
  N1s!0 = "run"
  N1s&run = runs
  N1s!1 = "lev"
  N1s&lev = lev_p
  N1s!2 = "lat"
  N1s&lat = lat
  N1s@long_name = "Concentration, r >= 1 nm"
  N1s@units = "# cm~S~-3~N~"  

  
 ; Calculate the zonal monthly average of particles r > 1 nm
  print("Averaging July N1 ...")
  T_month=dim_avg_Wrap(totsul10A(lev|:,lat|:,lon|:,time|time0701:time0731))
;  T_month=dim_avg_Wrap(pursul10A(lev|:,lat|:,lon|:,time|time0701:time0731))
;  T_month=dim_avg_Wrap(mixsul10A(lev|:,lat|:,lon|:,time|time0701:time0731))
  print("Zonally averaging July N1 ...")  
  N1s(lev|:,lat|:, run|1) = dim_avg(T_month(lev|:,lat|:,lon|:))
  
  print("Averaging January N1 ...")
  T_month=dim_avg_Wrap(totsul10A2(lev|:,lat|:,lon|:,time|time0101:time0131))
;  T_month=dim_avg_Wrap(pursul10A2(lev|:,lat|:,lon|:,time|time0101:time0131))
;  T_month=dim_avg_Wrap(mixsul10A2(lev|:,lat|:,lon|:,time|time0101:time0131))
  print("Zonally averaging January N1 ...")  
  N1s(lev|:,lat|:, run|0) = dim_avg(T_month(lev|:,lat|:,lon|:))
  

  ; Make a 4 plot panel ...
  print("")
  print("  Plotting ...")
 
  minP      = 0.0001        ; largest pressure to plot ...
  maxP      = 0.1        ; largest pressure to plot ...
  
  nPlots = 2
  plot = new((/ nPlots /), graphic)
  plotT = new((/ nPlots /), graphic)
  
  resP                      = True
  resP@gsnMaximize          = True
  resP@gsnPanelLabelBar     = True
;  resP@lbTitleOn            =  True                ; turn on title
;  resP@lbTitlePosition      = "Right"              ; title position
;  resP@lbTitleFontHeightF   = .015                ; make title smaller
;  resP@lbTitleDirection     = "Across"             ; title direction
  resP@gsnPanelXWhiteSpacePercent = 0.5
  resP@gsnPanelYWhiteSpacePercent = 0.5
  resP@lbLabelStride        = 1
  resP@lbLabelAngleF        = 60
  resP@tiMainString          = runname
        
  res                       = True
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@lbLabelBarOn          = False
  res@cnLinesOn             = False
  res@cnLinesOn             = False
  res@cnFillOn              = True                  ; turn on color fill
  res@cnInfoLabelOn         = False
  res@cnLineLabelsOn        = False
  res@gsnSpreadColors       = True                  ; use full range of colors
  res@cnFixFillBleed        = True
;  res@cnFillMode            = "CellFill"
;  res@lbTitlePosition       = "Right"              ; title position
;  res@lbTitleDirection      = "Across"             ; title direction
;  resP@lbTitleFontHeightF   = .015                ; make title smaller

  res@tmYLMode      = "Explicit" ; Define own tick mark labels.
;  res@tmYLLabelFont = 25         ; Change the font.
  res@tmYLValues    = (/0.1, 0.01, 0.001, 0.0001/)
  res@tmYLLabels    = (/"0.1", "0.01", "0.001", "0.0001"/)
  res@tmYLMinorOn   = True        ; No minor tick marks

  res@tmXBLabelFontHeightF  = 0.03   ; increase tickmark font size
  res@tmYRLabelFontHeightF  = 0.03   ; increase tickmark font size
  res@tmYLLabelFontHeightF  = 0.03   ; increase tickmark font size
  res@tmYLLabelAngleF      = 60
  
;  res@tmYRMode      = "Automatic" ; Define own tick mark labels.
;  res@tmYRLabelFont = 25         ; Change the font.
;  res@tmYRValues    = (/70.,80.,90.,100./)
;  res@tmYRLabels    = (/"70","80","90","100"/)
;  res@tmYRMinorOn   = True        ; No minor tick marks
    
  resO                      = True
  resO@gsnDraw              = False
  resO@gsnFrame             = False
  resO@lbLabelBarOn         = False
  resO@cnLineColor          = "Black"                   ; color of second contours
  resO@gsnRightString       = ""    
  resO@gsnLeftString        = ""
  resO@gsnCenterString      = ""
  resO@cnLevelSelectionMode = "ManualLevels" 
  resO@cnLevelSpacingF      = 10.0                   ; contour interval
  resO@cnMinLevelValF       = 110.0             ; min level
  resO@cnMaxLevelValF       = 160.0            ; max level
  resO@cnLineThicknessF     = 2.0                   ; color of second contours

  months = (/ "Dec", "January", "Feb", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov" /)

  ; Open a workstation
  wks = gsn_open_wks("pdf","images/totsulf125")
  gsn_define_colormap(wks,"gui_default")                ; choose colormap
  
  print("")
  print("    temperature overlay")
  zonal = Ts
  
  do i = 0, nPlots-1
    plotT(i) = gsn_csm_contour(wks, zonal(i, {minP:maxP}, :), resO)          ; draw second plot
;    plotT(i) = gsn_csm_contour(wks, zonal(i, :, :), resO)          ; draw second plot
  end do

  print("")
  print("    number density")
  zonal                     = N1s
  
;  res@lbTitleString        = zonal@units
;  res@lbLabelStride         = 3

  res@cnLevelSelectionMode  = "ManualLevels" 
;  print("      min :  " + min(zonal(:, {minP:maxP}, :)))
;  print("      max :  " + max(zonal(:, {minP:maxP}, :)))

; Contours by 10 up to 200:
;  res@cnLevelSpacingF       = 10.0                  ; contour interval
;  res@cnMinLevelValF        = 10.0                  ; contour interval
;  res@cnMaxLevelValF        = 200.0                 ; contour interval

;;  res@tiMainString          = runname

  res@cnLevelSelectionMode  = "ExplicitLevels"
  res@cnLevels = (/ 5, 10, 20, 30, 50, 75, 100, \
      200, 300, 500, 750, 1000, 2000, 3000, 5000 /);, 7500, 10000 /)
  res@lbLabelAngleF = 60
  
  res@tiMainString          = months(month(0))
  res@tiMainFontHeightF         = 0.04
  
  plot(0) = gsn_csm_pres_hgt(wks, zonal(0, {minP:maxP}, {-90.0:-30.0}), res)
;  plot(0) = gsn_csm_pres_hgt(wks, zonal(0, {minP:maxP}, {30.0:90.0}), res)
;  plot(0) = gsn_csm_pres_hgt(wks, zonal(0, :, {-90.0:-30.0}), res)
  overlay(plot(0), plotT(0))

  res@tiMainString          = months(month(1))
;  res@tiMainString          = runs(1)
;  plot(1) = gsn_csm_pres_hgt(wks, zonal(1, {minP:maxP}, {30.0:90.0}), res)
  plot(1) = gsn_csm_pres_hgt(wks, zonal(1, {minP:maxP}, {30.0:90.0}), res)
  overlay(plot(1), plotT(1))
    

;  print("")
;  print("    mass density")
;  zonal = MDs
;  zonal&lev = lev_p

;  res@lbTitleString        = zonal@units
;  print("      min :  " + min(zonal(:, {0.0:maxP}, :)))
;  print("      max :  " + max(zonal(:, {0.0:maxP}, :)))

;  res@cnLevelSelectionMode = "ExplicitLevels"
;  res@cnLevels = (/ .2, .4, .8, 1, 2, 4, 8, 10, 20, 40, 80, 100, 200 /)
;  res@lbLabelStride        = 1

;  res@tiMainString          = ""
;  do i = 0, nRuns-1
;    plot(2+i) = gsn_csm_pres_hgt(wks, zonal(i, {0.0:maxP}, :), res)
;  end do

  gsn_panel(wks, plot, (/ 1, 2 /), resP)

end


monave_za()
;monave_polar()
