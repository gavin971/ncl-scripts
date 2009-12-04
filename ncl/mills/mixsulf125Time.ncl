; Dust-Sulfates
;
; Daily average plots of aerosols larger than various radius thresholds:
;  1 nm (10A), 0.5 nm (05A), or all aerosols (00A)
;
; The routine analysisSulfLimit.ncl must be run first to generate the files
; in ~/analysis/<caseid> used as input.  Then the routine mixsulTimeAnalysis
; (below) is called first to produce the file mixsulf125Time.nc from these
; files.  This file is then used as input for the plotting routine,
; mixsulf125Time (below).

load "tests/common.ncl"
;load "tests/getHistoryPath.ncl"

undef("mixsulfTimeAnalysis")
procedure mixsulfTimeAnalysis(caseid:string,pattern:string)

local files, nFiles, historyFile

begin
  print("")
  print("mixsulfTimeAnalysis: caseid="+caseid)
  
; Determine the path to analysis files
  analysisPath = "analysis/" + caseid+ "/"

  ; Open the  history files
  analysisNames = systemfunc("ls analysis/" + caseid + "/*" + pattern + "*")
  fileNames = systemfunc("ls " + analysisPath + "| grep " + pattern)
    
  analysis = addfiles(analysisNames, "r")
  ListSetType(analysis, "cat")
  print("analysisNames:"+dimsizes(analysisNames))
   
 ; Compute the pressure levels ...
  time    = addfiles_GetVar(analysis, analysisNames, "time")
  nTime   = dimsizes(time)

  lev     = analysis[4]->lev
  nLevs   = dimsizes(lev)
  lat     = analysis[0]->lat
  nLats   = dimsizes(lat)
  lon     = analysis[0]->lon
  nLons   = dimsizes(lon)
  
  ; Determine an approximate dz/dp.
;  hyam 	= analysis[0]->hyam
;  hybm 	= analysis[0]->hybm
;  hyai 	= analysis[0]->hyai
;  hybi 	= analysis[0]->hybi
;  P0 	  = analysis[0]->P0
;  PS 	  = analysis[0]->PS
  
  print("reading date")
  date = addfiles_GetVar(analysis, analysisNames, "date")
  
  print("reading zonal mean data")
  T_zm  = addfiles_GetVar(analysis, analysisNames, "T_zm")
  pursulnd_zm  = addfiles_GetVar(analysis, analysisNames, "pursulnd_zm")
  mixsulnd_zm  = addfiles_GetVar(analysis, analysisNames, "mixsulnd_zm")
  dustnd_zm  = addfiles_GetVar(analysis, analysisNames, "dustnd_zm")

  mixsul05a_zm  = addfiles_GetVar(analysis, analysisNames, "mixsul05a_zm")
  mixsul10a_zm  = addfiles_GetVar(analysis, analysisNames, "mixsul10a_zm")
  mixsul13a_zm  = addfiles_GetVar(analysis, analysisNames, "mixsul13a_zm")
  mixsul30a_zm  = addfiles_GetVar(analysis, analysisNames, "mixsul30a_zm")

  pursul05a_zm  = addfiles_GetVar(analysis, analysisNames, "pursul05a_zm")
  pursul10a_zm  = addfiles_GetVar(analysis, analysisNames, "pursul10a_zm")
  pursul13a_zm  = addfiles_GetVar(analysis, analysisNames, "pursul13a_zm")
  pursul30a_zm  = addfiles_GetVar(analysis, analysisNames, "pursul30a_zm")

  dust05a_zm  = addfiles_GetVar(analysis, analysisNames, "dust05a_zm")
  dust10a_zm  = addfiles_GetVar(analysis, analysisNames, "dust10a_zm")
  dust13a_zm  = addfiles_GetVar(analysis, analysisNames, "dust13a_zm")
  dust30a_zm  = addfiles_GetVar(analysis, analysisNames, "dust30a_zm")
  
  print("reading 0.007 hPa data")
  T_007hPa  = addfiles_GetVar(analysis, analysisNames, "T_007hPa")
  pursulnd_007hPa  = addfiles_GetVar(analysis, analysisNames, "pursulnd_007hPa")
  mixsulnd_007hPa  = addfiles_GetVar(analysis, analysisNames, "mixsulnd_007hPa")
  dustnd_007hPa  = addfiles_GetVar(analysis, analysisNames, "dustnd_007hPa")

  mixsul05a_007hPa  = addfiles_GetVar(analysis, analysisNames, "mixsul05a_007hPa")
  mixsul10a_007hPa  = addfiles_GetVar(analysis, analysisNames, "mixsul10a_007hPa")
  mixsul13a_007hPa  = addfiles_GetVar(analysis, analysisNames, "mixsul13a_007hPa")
  mixsul30a_007hPa  = addfiles_GetVar(analysis, analysisNames, "mixsul30a_007hPa")

  pursul05a_007hPa  = addfiles_GetVar(analysis, analysisNames, "pursul05a_007hPa")
  pursul10a_007hPa  = addfiles_GetVar(analysis, analysisNames, "pursul10a_007hPa")
  pursul13a_007hPa  = addfiles_GetVar(analysis, analysisNames, "pursul13a_007hPa")
  pursul30a_007hPa  = addfiles_GetVar(analysis, analysisNames, "pursul30a_007hPa")

  dust05a_007hPa  = addfiles_GetVar(analysis, analysisNames, "dust05a_007hPa")
  dust10a_007hPa  = addfiles_GetVar(analysis, analysisNames, "dust10a_007hPa")
  dust13a_007hPa  = addfiles_GetVar(analysis, analysisNames, "dust13a_007hPa")
  dust30a_007hPa  = addfiles_GetVar(analysis, analysisNames, "dust30a_007hPa")
  
  lev_p   = lev
  lev_p@units = "mb"
  lev_p@long_name = "Pressure (mb)" 
  
 ; Calculate the zonal daily average of temperature
  T_zm&lev = lev_p  
  pursulnd_zm&lev = lev_p
  mixsulnd_zm&lev = lev_p
  dustnd_zm&lev = lev_p
  
  pursul05a_zm&lev = lev_p  
  pursul10a_zm&lev = lev_p  
  pursul13a_zm&lev = lev_p  
  pursul30a_zm&lev = lev_p
  
  dust05a_zm&lev = lev_p  
  dust10a_zm&lev = lev_p  
  dust13a_zm&lev = lev_p  
  dust30a_zm&lev = lev_p
  
  mixsul05a_zm&lev = lev_p  
  mixsul10a_zm&lev = lev_p  
  mixsul13a_zm&lev = lev_p  
  mixsul30a_zm&lev = lev_p
  
  if (.not.(isfilepresent("analysis/"+caseid))) then
    a = systemfunc("mkdir -p analysis/"+caseid)
  end if
  
  outfile="analysis/"+caseid+"/mixsulf125Time."+caseid+".nc"

  print("saving file "+outfile)
  fileout=addfile(outfile,"c")
  
  fileout->date=date
  fileout->T_zm=T_zm
  fileout->T_007hPa=T_007hPa
  
  fileout->mixsulnd_zm=mixsulnd_zm
  fileout->mixsul05a_zm=mixsul05a_zm
  fileout->mixsul10a_zm=mixsul10a_zm
  fileout->mixsul13a_zm=mixsul13a_zm
  fileout->mixsul30a_zm=mixsul30a_zm
  fileout->mixsulnd_007hPa=mixsulnd_007hPa
  fileout->mixsul05a_007hPa=mixsul05a_007hPa
  fileout->mixsul10a_007hPa=mixsul10a_007hPa
  fileout->mixsul13a_007hPa=mixsul13a_007hPa
  fileout->mixsul30a_007hPa=mixsul30a_007hPa
  
  fileout->pursulnd_zm=pursulnd_zm
  fileout->pursul05a_zm=pursul05a_zm
  fileout->pursul10a_zm=pursul10a_zm
  fileout->pursul13a_zm=pursul13a_zm
  fileout->pursul30a_zm=pursul30a_zm
  fileout->pursulnd_007hPa=pursulnd_007hPa
  fileout->pursul05a_007hPa=pursul05a_007hPa
  fileout->pursul10a_007hPa=pursul10a_007hPa
  fileout->pursul13a_007hPa=pursul13a_007hPa
  fileout->pursul30a_007hPa=pursul30a_007hPa
  
  fileout->dustnd_zm=dustnd_zm
  fileout->dust05a_zm=dust05a_zm
  fileout->dust10a_zm=dust10a_zm
  fileout->dust13a_zm=dust13a_zm
  fileout->dust30a_zm=dust30a_zm
  fileout->dustnd_007hPa=dustnd_007hPa
  fileout->dust05a_007hPa=dust05a_007hPa
  fileout->dust10a_007hPa=dust10a_007hPa
  fileout->dust13a_007hPa=dust13a_007hPa
  fileout->dust30a_007hPa=dust30a_007hPa

end

undef("mixsulfTime")
procedure mixsulfTime(caseid:string,targetLat:float,year:string,sizeLimit:string)

; inputs:
;  caseid, e.g. "125kelvH2Oe"
;  targetLat, e.g. 84 for 84°N and 84°S
;  year = can be scalar string year for NH only, or
;         array of (/yearNH,yearSH1,yearSH2/) for 2-hemisphere plots

local files, nFiles, inFile

begin
  
  print("")
  print("caseid="+caseid)
  print("  Plot daily average mixed sulfate densities ...")
  
  ; Open the  history files
  
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
    dust_zm   = filein->dust10a_zm
  else
    if (sizeLimit.eq."05A") then
      mixsul_zm = filein->mixsul05a_zm
      pursul_zm = filein->pursul05a_zm
      dust_zm   = filein->dust05a_zm
    else
      if (sizeLimit.eq."00A") then  
        mixsul_zm = filein->mixsulnd_zm
        pursul_zm = filein->pursulnd_zm
        dust_zm   = filein->dustnd_zm
      else
        print("FATAL: sizeLimit not recognized="+sizeLimit)
        return
      end if
    end if
  end if
        
  nYear = dimsizes(year)
  print("nYear="+nYear)
  
  time0401=ind(date.eq.year(0)+"0401")
  time1031=ind(date.eq.year(0)+"1031")
  
  if (ismissing(time1031)) then
    time1031=nTime-1
  end if
  
  print("time0401="+time0401+" time="+time(time0401))
  print("time1031="+time1031+" time="+time(time1031))
 
  if (nYear.gt.1) then
  
    time1001=ind(date.eq.year(1)+"1001")
    time0430=ind(date.eq.year(2)+"0430")

    print("time1001="+time1001+" time="+time(time1001))
    print("time0430="+time0430+" time="+time(time0430))
  
  else ; fill in some values for tick generation code below
  
    time1001=time0401
    time0430=time1031

  end if
  
  ; Make a 3- or 4-plot panel ...
  print("")
  print("  Plotting ...")
 
  minP      = 0.0001        ; largest pressure to plot ...
  maxP      = 0.1        ; largest pressure to plot ...
  
  nPlots = 3
  if (nYear.gt.1) then
    nPlots = 4
  end if
  
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
  
  if (.not.(isfilepresent("images/"+caseid))) then
    a = systemfunc("mkdir -p images/"+caseid)
  end if
  
  wks = gsn_open_wks("pdf","images/"+caseid+"/mixsulf125Time."+caseid+"."+sizeLimit+"."+doubletoint(lat({targetLat})))
  gsn_define_colormap(wks,"gui_default")                ; choose colormap
  
  print("")
  print("    temperature overlay")
  
 ; Create temperature overlay contours
  plotT(0) = gsn_csm_contour(wks, T_zm({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat}), resO)
  
  if (nPlots.eq.3) then
  
    plotT(1) = gsn_csm_contour(wks, T_zm({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat}), resO)  
    plotT(2) = gsn_csm_contour(wks, T_zm({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat}), resO)  
    
  else
  
    plotT(2) = gsn_csm_contour(wks, T_zm({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat}), resO)

    plotT(1) = gsn_csm_contour(wks, T_zm({lev|minP:maxP} , time|time1001:time0430, {lat|-targetLat}), resO)
    plotT(3) = gsn_csm_contour(wks, T_zm({lev|minP:maxP} , time|time1001:time0430, {lat|-targetLat}), resO)
    
  end if
  
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
;  res@cnLevels = (/ 5, 10, 20, 30, 50, 75, 100, \
;      200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000, 7500, 10000 /)
  res@cnLevels = (/ 1, 2, 4, 8, 10, 20, 40, 80, 100, 200, 400, 600, 800, 1000, 2000, 4000, 8000, 1e4 /)
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
  
  plot(0) = gsn_csm_pres_hgt(wks, mixsul_zm({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat}), res)
  overlay(plot(0), plotT(0))

  if (nPlots.eq.3) then
  
    plot(1) = gsn_csm_pres_hgt(wks, pursul_zm({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat}), res)
    overlay(plot(1), plotT(1))
  
    plot(2) = gsn_csm_pres_hgt(wks, dust_zm({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat}), res)
    overlay(plot(2), plotT(2))
    
    gsn_panel(wks, plot, (/ 1, 3 /), resP)
  
  else
    plot(1) = gsn_csm_pres_hgt(wks, mixsul_zm({lev|minP:maxP} , time|time1001:time0430, {lat|-targetLat}), resSH)
    overlay(plot(1), plotT(1))

    plot(2) = gsn_csm_pres_hgt(wks, pursul_zm({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat}), res)
    overlay(plot(2), plotT(2))

    plot(3) = gsn_csm_pres_hgt(wks, pursul_zm({lev|minP:maxP} , time|time1001:time0430, {lat|-targetLat}), resSH)
    overlay(plot(3), plotT(3))

    gsn_panel(wks, plot, (/ 2, 2 /), resP)
  end if
  
  print("NH mixsulf max: "+max(mixsul_zm({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat})))
  print("NH pursulf max: "+max(pursul_zm({lev|minP:maxP} , time|time0401:time1031, {lat|targetLat})))

end

caseid="125b36Pb02"
year="1996"
lat=86

;mixsulfTimeAnalysis(caseid,"h3")

;mixsulfTime(caseid,70)
mixsulfTime(caseid,lat,year,"00A")
mixsulfTime(caseid,lat,year,"05A")
mixsulfTime(caseid,lat,year,"10A")
