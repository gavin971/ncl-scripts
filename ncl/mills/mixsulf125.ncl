; Dust-Sulfates
;
; Monthly average plots

load "tests/common.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("monave_za")
procedure monave_za(caseid:string, aerType:string, year:string, sizeLimit:string)

;  year = can be scalar string year for NH only, or
;         array of (/yearNH,yearSH1,yearSH2/) for 2-hemisphere plots

local files, nFiles, historyFile

begin
  
  print("")
  print("caseid="+caseid)
  print("  Plot monthly average mixed sulfate densities ...")

  ; Open the analysis file
  infile="analysis/"+caseid+"/mixsulf125Time."+caseid+".nc"    
  filein=addfile(infile,"r")

  time    = filein->time
  nTime   = dimsizes(time)

  lev     = filein->lev
  nLevs   = dimsizes(lev)
  lat     = filein->lat
  nLats   = dimsizes(lat)
  
  date    = filein->date
  T_zm    = filein->T_zm

; Set for various size thresholds

  if (sizeLimit.eq."10A") then
    mixsul_zm = filein->mixsul10a_zm
    pursul_zm = filein->pursul10a_zm
  else
    if (sizeLimit.eq."05A") then
      mixsul_zm = filein->mixsul05a_zm
      pursul_zm = filein->pursul05a_zm
    else
      if (sizeLimit.eq."00A") then  
        mixsul_zm = filein->mixsulnd_zm
        pursul_zm = filein->pursulnd_zm
      else
        print("FATAL: sizeLimit not recognized="+sizeLimit)
        return
      end if
    end if
  end if
        
  nYear = dimsizes(year)
  print("nYear="+nYear)
  
  yearNH = year(0)
  startdate=yearNH+"0701"
  enddate=yearNH+"0731"
  print(startdate)
  print(enddate)
  time1s=ind(date.eq.stringtoint(startdate))
  time1e=ind(date.eq.stringtoint(enddate))
  print("time1s: "+startdate+" | "+time1s)
  print("time1e: "+enddate+" | "+time1e)
  
  if (nYear .gt. 1) then
  
    yearSH1= year(1)
    yearSH2= year(2)      
    startdate=yearSH2+"0101"
    enddate=yearSH2+"0131"
    print(startdate)
    print(enddate)
    time2s=ind(date.eq.stringtoint(startdate))
    time2e=ind(date.eq.stringtoint(enddate))
    print("time2s: "+startdate+" | "+time2s)
    print("time2e: "+enddate+" | "+time2e)
    
    month = (/ 1, 7 /)
    runs = (/ "January", "July" /)

    if (aerType .eq. "mix") then
      sulfate1_zm  = mixsul_zm
      sulfate2_zm  = mixsul_zm
    end if

    if (aerType .eq. "pur") then
      sulfate1_zm  = pursul_zm
      sulfate2_zm  = pursul_zm
    end if

    if (aerType .eq. "tot") then
      sulfate1_zm  = mixsul_zm ; get dimension wrappers
      sulfate1_zm  = sulfate1_zm + pursul_zm
      sulfate2_zm  = sulfate1_zm
    end if

    lat21=-90.0
    lat22=-30.0
    
  else
  
    month = (/ 7, 7 /)
    runs = (/ "Mixed Sulfates", "Pure Sulfates" /)
    
    time2s=time1s
    time2e=time1e
    
    sulfate1_zm  = mixsul_zm
    sulfate2_zm  = pursul_zm
    
    lat21=30.0
    lat22=90.0

  end if
  
  lev_p   = lev
  lev_p@units = "mb"
  lev_p@long_name = "Pressure (mb)"
  
  nMonths = dimsizes(month)
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
  print("Averaging monthly temperature ...")
  Ts(lev|:,lat|:, run|0) = dim_avg(T_zm(lev|:,lat|:,time|time1s:time1e))
  Ts(lev|:,lat|:, run|1) = dim_avg(T_zm(lev|:,lat|:,time|time2s:time2e))
  
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
  print("Averaging monthly sulfate  ...")
  N1s(lev|:,lat|:, run|0)=dim_avg_Wrap(sulfate1_zm(lev|:,lat|:,time|time1s:time1e))
  N1s(lev|:,lat|:, run|1)=dim_avg_Wrap(sulfate2_zm(lev|:,lat|:,time|time2s:time2e))
  
  ; Make a 2-plot panel ...
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
  resP@tiMainString          = caseid
        
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

  ; Open a workstation
  wks = gsn_open_wks("pdf","images/"+caseid+"/"+aerType+"sulf125."+caseid+"."+sizeLimit)
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


  res@cnLevelSelectionMode  = "ExplicitLevels"
  res@cnLevels = (/ 5, 10, 20, 30, 50, 75, 100, \
      200, 300, 400, 500, 600, 700, 800/)
  res@lbLabelAngleF = 60
  
  res@tiMainString          = runs(0)
  res@tiMainFontHeightF         = 0.04
  
  plot(0) = gsn_csm_pres_hgt(wks, zonal(0, {minP:maxP}, {lat21:lat22}), res)
;  plot(0) = gsn_csm_pres_hgt(wks, zonal(0, {minP:maxP}, {30.0:90.0}), res)
;  plot(0) = gsn_csm_pres_hgt(wks, zonal(0, :, {-90.0:-30.0}), res)
  overlay(plot(0), plotT(0))

  res@tiMainString          = runs(1)
;  res@tiMainString          = runs(1)
  plot(1) = gsn_csm_pres_hgt(wks, zonal(1, {minP:maxP}, {30.0:90.0}), res)
  overlay(plot(1), plotT(1))
    
  print("      min :  " + min(zonal(1, {minP:maxP}, {60.0:90.0})))
  print("      max :  " + max(zonal(1, {minP:maxP}, {60.0:90.0})))

  gsn_panel(wks, plot, (/ 1, 2 /), resP)

end ; procedure monave_za(caseid:string, aerType:string, year:string, sizeLimit:string)

caseid="125b36eg"

;aerType = "tot"
aerType = "mix"
;aerType = "pur"

monave_za(caseid,aerType,"1995","00A")
monave_za(caseid,aerType,"1995","05A")
monave_za(caseid,aerType,"1995","10A")
;monave_polar(caseid)
