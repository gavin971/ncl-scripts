; Dust-Sulfates
;
; Daily average plots

load "tests/common.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("mixsulfTimeAnalysis")
procedure mixsulfTimeAnalysis(caseid:string)

local files, nFiles, historyFile

begin
  print("")
  print("caseid="+caseid)
  print("  Plot daily average mixed sulfate densities ...")
  
  ; Open the  history files
  
  historyPattern="waccm/dust-sulfur/run/"+caseid+"/alldaily/all.*.cam2.h3.*"    
  
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  print("historyNames:"+dimsizes(historyNames))
;  print(historyNames)
   
 ; Compute the pressure levels ...
  time    = addfiles_GetVar(history, historyNames, "time")
  nTime   = dimsizes(time)

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
  
  print("reading date")
  date = addfiles_GetVar(history, historyNames, "date")
  print("reading T")
  T    = addfiles_GetVar(history, historyNames, "T")
  
  print("reading mixsul10A")
  mixsul10A  = addfiles_GetVar(history,  historyNames,  "MIXSUL10A")
  print("reading pursul10A")
  pursul10A  = addfiles_GetVar(history,  historyNames,  "PURSUL10A")
;  totsul10A  = mixsul10A ; get dim wrappers
;  totsul10A  = totsul10A + pursul10A
  
  
  lev_p   = lev
  lev_p@units = "mb"
  lev_p@long_name = "Pressure (mb)" 
  
 ; Calculate the zonal daily average of temperature
  print("Averaging temperature ...")
  T_za=dim_avg_Wrap(T(lev|:,time|:,lat|:,lon|:))
  T_za&lev = lev_p
  
  print("Averaging pursulf ...")
  pursulf_za=dim_avg_Wrap(pursul10A(lev|:,time|:,lat|:,lon|:))
  pursulf_za&lev = lev_p
  
  print("Averaging mixsulf ...")
  mixsulf_za=dim_avg_Wrap(mixsul10A(lev|:,time|:,lat|:,lon|:))  
  mixsulf_za&lev = lev_p
  
  outfile="analysis/"+caseid+"/mixsulf125Time.nc"

  print("saving file "+outfile)
  fileout=addfile(outfile,"c")
  
  print("saving date")
  fileout->date=date
  print("saving T_za")
  fileout->T_za=T_za
  print("saving mixsulf_za")
  fileout->mixsulf_za=mixsulf_za
  print("saving pursulf_za")
  fileout->pursulf_za=pursulf_za

end

undef("mixsulfTime")
procedure mixsulfTime(caseid:string,targetLat:float)

local files, nFiles, historyFile

begin
  
  print("")
  print("caseid="+caseid)
  print("  Plot daily average mixed sulfate densities ...")
  
  ; Open the  history files
  
  infile="analysis/"+caseid+"/mixsulf125Time.nc"    
  filein=addfile(infile,"r")

  time    = filein->time
  nTime   = dimsizes(time)

  lev     = filein->lev
  nLevs   = dimsizes(lev)
  lat     = filein->lat
  nLats   = dimsizes(lat)
  
  date    = filein->date
  T_za    = filein->T_za
  mixsulf_za = filein->mixsulf_za
  pursulf_za = filein->pursulf_za
  
  yearNH="1997"
  yearSH1="1997"
  yearSH2="1998"
  
  if (caseid .eq. "125smin400" .or caseid .eq. "125tabaz22") then
    yearNH="1995"
    yearSH1="1995"
    yearSH2="1996"
  end if
  
  if (caseid .eq. "125sminPtbo2") then
    yearNH="1996"
    yearSH1="1995"
    yearSH2="1996"
  end if
  
  time0401=ind(date.eq.yearNH+"0401")
  time1031=ind(date.eq.yearNH+"1031")
  
  print("time0401="+time0401+" time="+time(time0401))
  print("time1031="+time1031+" time="+time(time1031))
 
  time1001=ind(date.eq.yearSH1+"1001")
  time0430=ind(date.eq.yearSH2+"0430")
  
  print("time1001="+time1001+" time="+time(time1001))
  print("time0430="+time0430+" time="+time(time0430))

  ; Make a 4 plot panel ...
  print("")
  print("  Plotting ...")
 
  minP      = 0.0001        ; largest pressure to plot ...
  maxP      = 0.1        ; largest pressure to plot ...
  
  nPlots = 4
  plot = new((/ nPlots /), graphic)
  plotT = new((/ nPlots /), graphic)
  
  resP                      = True
  resP@gsnMaximize          = True
  resP@gsnBoxMargin          = 0.0
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
  
;  res@tmYRMode      = "Automatic" ; Define own tick mark labels.
;  res@tmYRLabelFont = 25         ; Change the font.
;  res@tmYRValues    = (/70.,80.,90.,100./)
;  res@tmYRLabels    = (/"70","80","90","100"/)
;  res@tmYRMinorOn   = True        ; No minor tick marks
  res@tmXBLabelFontHeightF  = 0.05   ; increase tickmark font size
  res@tmYRLabelFontHeightF  = 0.04   ; increase tickmark font size
  res@tmYLLabelFontHeightF  = 0.03   ; increase tickmark font size
  res@tmYLLabelAngleF      = 60
  res@gsnMaximize           = True
  res@gsnBoxMargin          = 0.0
  res@tiXAxisFontHeightF = 0.027
  res@tiYAxisFontHeightF = 0.027
  
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
  resO@gsnMaximize          = True
  resO@gsnBoxMargin          = 0.0

  months = (/ "Dec", "January", "Feb", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov" /)

  ; Open a workstation
  wks = gsn_open_wks("pdf","images/mixsulf125Time."+caseid+"."+targetLat)
  gsn_define_colormap(wks,"gui_default")                ; choose colormap
  
  print("")
  print("    temperature overlay")
  
 ; Create temperature overlay contours
  plotT(0) = gsn_csm_contour(wks, T_za({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat}), resO)
  plotT(2) = gsn_csm_contour(wks, T_za({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat}), resO)
  
  plotT(1) = gsn_csm_contour(wks, T_za({lev|minP:maxP} , time|time1001:time0430, {lat|-targetLat}), resO)
  plotT(3) = gsn_csm_contour(wks, T_za({lev|minP:maxP} , time|time1001:time0430, {lat|-targetLat}), resO)

  print("")
  print("    number density")
   
;  res@lbTitleString        = zonal@units
;  res@lbLabelStride         = 3

;  res@cnLevelSelectionMode  = "ManualLevels" 
;  print("      min :  " + min(zonal(:, {minP:maxP}, :)))
;  print("      max :  " + max(zonal(:, {minP:maxP}, :)))
;  res@cnLevelSpacingF       = 5.0                  ; contour interval
;  res@cnMinLevelValF        = 5.0                  ; contour interval
;  res@cnMaxLevelValF        = 100.0                 ; contour interval  

  res@cnLevelSelectionMode  = "ExplicitLevels"
  res@cnLevels = (/ 5, 10, 20, 30, 50, 75, 100, \
      200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000, 7500, 10000 /)
  res@lbLabelAngleF = 60
  
;===================================================================
; Do special labels on X axis. There will be 6 tickmarks per each
; month, with a month label appearing in the middle.  Here, we are
; treating minor tickmarks like major tickmarks, and vice versa, since
; only the major tickmarks can actually be labeled.
;
  nmonths =8
  
  monthsNH = ispan(4,11,1)
  
  monthsSH = new(nmonths,integer)
  monthsSH(0:2) = ispan(10,12,1)
  monthsSH(3:7) = ispan(1,5,1)

  day_in_monNH  = new(nmonths,float)           ; Number of days in each month
  day_in_monSH  = new(nmonths,float)           ; Number of days in each month
  do i=0,nmonths-1
    day_in_monNH(i) = days_in_month(1997,monthsNH(i))*1.
    day_in_monSH(i) = days_in_month(1997,monthsSH(i))*1.
  end do

  imonth      = 0
  days_accumNH  = 15.+stringtofloat(time(time0401))
  days_accumSH  = 15.+stringtofloat(time(time1001))

  tk_per_mon  = 4                    ; Tickmarks per month.
  nticks      = nmonths * tk_per_mon
  major_ticksNH = new(nticks,float)    ; Array to hold major tickmark values.
  major_ticksSH = new(nticks,float)    ; Array to hold major tickmark values.
  ii5         = ispan(1,tk_per_mon-1,1)*1.
  do i=0,nticks-1,tk_per_mon
    major_ticksNH(i) = days_accumNH
    major_ticksSH(i) = days_accumSH
    major_ticksNH(i+1:i+tk_per_mon-1) = major_ticksNH(i) + ii5*(day_in_monNH(imonth)/tk_per_mon)
    major_ticksSH(i+1:i+tk_per_mon-1) = major_ticksSH(i) + ii5*(day_in_monSH(imonth)/tk_per_mon)
    days_accumNH     = days_accumNH + day_in_monNH(imonth)
    days_accumSH     = days_accumSH + day_in_monSH(imonth)
    imonth = imonth + 1
  end do
  res@tmXBMinorValues    = major_ticksNH(tk_per_mon/2:nticks-1:tk_per_mon)
  res@tmXBValues         = major_ticksNH
  res@tmXBLabels         = new(nticks,string)
  res@tmXBLabels(:)     = ""

  print("nticks="+nticks)
  monthnames=(/"J","F","M","A","M","J","J","A","S","O","N","D"/)
  print("monthnames:"+dimsizes(monthnames))
  print("labels:"+dimsizes(res@tmXBLabels(::tk_per_mon)))
  print("nticks/tk_per_mon="+nticks/tk_per_mon)
  
  res@tmXBMode           = "Explicit"         ; Define own tick mark labels.

  res@tmXBMinorOutwardLengthF = 0.02    ; Increase length and size of
  res@tmXBMinorLengthF        = 0.01    ; minor tickmarks
  res@tmXBMinorThicknessF     = 2.0
  res@tmXBMajorThicknessF     = 2.0
  res@tmXBMajorOutwardLengthF = 0.01
  res@tmXBMajorLengthF        = 0.01
  
 ; Label mid tickmark with the month label.
  res@tmXBLabels(::tk_per_mon)  = monthnames(monthsNH-1)
  
  resSH=res
  resSH@tmXBMinorValues    = major_ticksSH(tk_per_mon/2:nticks-1:tk_per_mon)
  resSH@tmXBValues         = major_ticksSH
  resSH@tmXBLabels         = new(nticks,string)
  resSH@tmXBLabels(:)     = ""
  resSH@tmXBLabels(::tk_per_mon)  = monthnames(monthsSH-1)
  
;===================================================================  
  
  plot(0) = gsn_csm_pres_hgt(wks, mixsulf_za({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat}), res)
  overlay(plot(0), plotT(0))

  plot(1) = gsn_csm_pres_hgt(wks, mixsulf_za({lev|minP:maxP} , time|time1001:time0430, {lat|-targetLat}), resSH)
  overlay(plot(1), plotT(1))
    
  plot(2) = gsn_csm_pres_hgt(wks, pursulf_za({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat}), res)
  overlay(plot(2), plotT(2))

  plot(3) = gsn_csm_pres_hgt(wks, pursulf_za({lev|minP:maxP} , time|time1001:time0430, {lat|-targetLat}), resSH)
  overlay(plot(3), plotT(3))


  gsn_panel(wks, plot, (/ 2, 2 /), resP)
  
  print("NH mixsulf max: "+max(mixsulf_za({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat})))
  print("NH pursulf max: "+max(pursulf_za({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat})))

end

caseid="125tabaz22"
;caseid="125smin400"
;caseid="125sminPtbo2"
;mixsulfTimeAnalysis(caseid)

;mixsulfTime(caseid,70)
mixsulfTime(caseid,84)
