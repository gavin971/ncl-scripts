; Dust-Sulfates
;
; Monthly average plots

load "tests/common.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("monave_za")
procedure monave_za(caseid1:string, caseid2:string)

local files, nFiles, historyFile

begin
  
  print("")
  caseidarr=(/caseid1,caseid2/)
  print("caseid1="+caseidarr(0))
  print("caseid2="+caseidarr(1))
  print("  Plot July average sulfate densities ...")

  ; Make a 4 plot panel ...
  print("")
  print("  Plotting ...")
 
  minP      = 0.0001        ; largest pressure to plot ...
  maxP      = 0.1        ; largest pressure to plot ...
  
  nPlots = 6
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
  resP@tiMainString          = caseid1+" "+caseid2
        
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
  wks = gsn_open_wks("pdf","images/sulf125."+caseid1+"."+caseid2)
  gsn_define_colormap(wks,"gui_default")                ; choose colormap

  
  ; Open the  history files
  
  do icase = 0,1
  
    caseid=caseidarr(icase)
    print(" case: "+caseid)

    if (caseid .eq. "125expuT15") then     
      yearNH="1997"
    end if

    if (caseid .eq. "125smin400") then
      yearNH="1995"
    end if

    if (caseid .eq. "125sminPtbo2") then
      yearNH="1996"
      yearSH1="1995"
      yearSH2="1996"
    end if

    if (caseid .eq. "125smin") then
      yearNH="1998"
    end if    

    historyPattern="waccm/dust-sulfur/run/"+caseid+"/alldaily/all.*.cam2.h3."+yearNH+"-0[67]*"
    historyNames = systemfunc("ls " + historyPattern)
    history = addfiles(historyNames, "r")
    ListSetType(history, "cat")
    print("historyNames:"+dimsizes(historyNames))
    print(historyNames)

   ; Compute the pressure levels ...
    time     = addfiles_GetVar(history, historyNames, "time")

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
    date = addfiles_GetVar(history, historyNames, "date")

    mixsul10A  = addfiles_GetVar(history,  historyNames,  "MIXSUL10A")
    pursul10A  = addfiles_GetVar(history,  historyNames,  "PURSUL10A")
    delete(history)
    delete(historyNames)
    
    totsul10A  = mixsul10A ; get dim wrappers
    totsul10A  = totsul10A + pursul10A

    print(date)
    startdate=yearNH+"0701"
    enddate=yearNH+"0731"
    print(startdate)
    print(enddate)
    time0701=ind(date.eq.stringtoint(startdate))
    time0731=ind(date.eq.stringtoint(enddate))
    print("time0701: "+startdate+" | "+time0701)
    print("time0731: "+enddate+" | "+time0731)

    lev_p   = lev
    lev_p@units = "mb"
    lev_p@long_name = "Pressure (mb)"

    runs = (/ "Total Sulfates","Pure Sulfates", "Mixed Sulfates" /)
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
    Ts(lev|:,lat|:, run|0) = dim_avg(T_month(lev|:,lat|:,lon|:))
    Ts(lev|:,lat|:, run|1) = Ts(lev|:,lat|:, run|0)
    Ts(lev|:,lat|:, run|2) = Ts(lev|:,lat|:, run|0)

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
    print("Averaging totsul10A July N1 ...")
    T_month=dim_avg_Wrap(totsul10A(lev|:,lat|:,lon|:,time|time0701:time0731))
    print("Zonally averaging totsul10A July N1 ...")  
    N1s(lev|:,lat|:, run|0) = dim_avg(T_month(lev|:,lat|:,lon|:))
    delete(totsul10A)

    print("Averaging pursul10A July N1 ...")
    T_month=dim_avg_Wrap(pursul10A(lev|:,lat|:,lon|:,time|time0701:time0731))
    N1s(lev|:,lat|:, run|1) = dim_avg(T_month(lev|:,lat|:,lon|:))
    delete(pursul10A)

    print("Averaging mixsul10A July N1 ...")
    T_month=dim_avg_Wrap(mixsul10A(lev|:,lat|:,lon|:,time|time0701:time0731))
    N1s(lev|:,lat|:, run|2) = dim_avg(T_month(lev|:,lat|:,lon|:))
    delete(mixsul10A)
    delete(time)
    delete(date)
    delete(T)

    print("")
    res@cnLevelSelectionMode  = "ManualLevels" 

    res@cnLevelSelectionMode  = "ExplicitLevels"
    res@cnLevels = (/ 5, 10, 20, 30, 50, 75, 100, \
        200, 300, 500, 750, 1000, 1500, 2000 , 3000, 5000, 7500, 10000 /)
    res@lbLabelAngleF = 60

    res@tiMainFontHeightF         = 0.04

    do i = 0, nRuns-1
      nplot=icase*nRuns+i
      res@tiMainString = runs(i)
      print("    nplot:"+nplot)
      plotT(nplot) = gsn_csm_contour(wks, Ts(i, {minP:maxP}, :), resO)          ; draw second plot
      plot(nplot) = gsn_csm_pres_hgt(wks, N1s(i, {minP:maxP}, {30.0:90.0}), res)
      overlay(plot(nplot), plotT(nplot))
    end do

  end do ; icase

  gsn_panel(wks, plot, (/ 2, 3 /), resP)

end

caseid1="125smin400"
caseid2="125sminPtbo2"


monave_za(caseid1,caseid2)
;monave_polar(caseid)
