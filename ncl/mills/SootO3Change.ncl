load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("plotO3Diff")
procedure plotO3Diff(O3_col_x[*][*]:float, O3_col_c[*][*]:float) 
local O3_col_diff,wks,res,lev
begin 

  lat=O3_col_x&lat
  
  xdim=dimsizes(O3_col_x)
  cdim=dimsizes(O3_col_c)
  print("xdim="+xdim+" cdim="+cdim)
  ntimes=min((/xdim(0),cdim(0)/))
  
  O3_col_diff=O3_col_x(:ntimes-1,:)-O3_col_c(:ntimes-1,:)
  O3_col_diff!0 = "time"
  O3_col_diff!1 = "lat"
  O3_col_diff&lat = lat

  O3_col_dpct=100.*O3_col_diff(:ntimes-1,:)/O3_col_c(:ntimes-1,:)
  O3_col_dpct!0 = "time"
  O3_col_dpct!1 = "lat"
  O3_col_dpct&lat = lat
  
  wks = gsn_open_wks("pdf", "O3ColumnChange")
  gsn_define_colormap(wks,"BlRe")
  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
;  res@lbLabelBarOn        = False         ; turn off color labelbar
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.55
  res@gsnLeftStringFontHeightF = 24
  res@gsnRightStringFontHeightF = 24
  res@tiYAxisString       = "Latitude"
  res@cnLineLabelPlacementMode= "Computed"
  res@cnLineLabelDensityF = 5
  res@tmXBLabelFontHeightF =0.015
  res@tmYLLabelFontHeightF =0.03
 
;===================================================================
; Do special labels on X axis. There will be 6 tickmarks per each
; month, with a month label appearing in the middle.  Here, we are
; treating minor tickmarks like major tickmarks, and vice versa, since
; only the major tickmarks can actually be labeled.
;
  year=2025
  yearlast=2035
  nyears=yearlast-year
  nmonths =8+12*nyears
  
  months  = new(nmonths,integer)
  years  = new(nmonths,integer)
  months(0:7) = ispan(5,12,1)
  years(0:7) = year
  do i=8,nmonths-1,12
    year=year+1
    years(i:i+11)=year
    months(i:i+11) =ispan(1,12,1) 
  end do

  day_in_mon  = new(nmonths,float)           ; Number of days in each month
  do i=0,nmonths-1
    day_in_mon(i) = days_in_month(years(i),months(i))*1.
  end do

  imonth      = 0
  days_accum  = 0.
  if (nmonths .lt. 60) then
    tk_per_mon  =4                                     ; Tickmarks per month.
    nticks       = nmonths * tk_per_mon
    major_ticks = new(nticks,float)    ; Array to hold major tickmark values.
    ii5         = ispan(1,tk_per_mon-1,1)*1.
    do i=0,nticks-1,tk_per_mon
      major_ticks(i) = days_accum
      major_ticks(i+1:i+tk_per_mon-1) = major_ticks(i) + ii5*(day_in_mon(imonth)/tk_per_mon)
      days_accum     = days_accum + day_in_mon(imonth)
      imonth = imonth + 1
    end do
    res@tmXBMinorValues    = major_ticks(tk_per_mon/2:nticks-1:tk_per_mon)
    res@tmXBValues         = major_ticks
    res@tmXBLabels         = new(nticks,string)
    res@tmXBLabels(:)     = ""

    print("nticks="+nticks)
    monthnames=(/"J","F","M","A","M","J","J","A","S","O","N","D"/)
    print("monthnames:"+dimsizes(monthnames))
    print("labels:"+dimsizes(res@tmXBLabels(::tk_per_mon)))
    print("nticks/tk_per_mon="+nticks/tk_per_mon)

   ; Label mid tickmark with the month label.
    res@tmXBLabels(::tk_per_mon)  = monthnames(months-1)
    
  else
    tk_per_year=12
    nticks=nmonths
    days_accum  = 0.  
    major_ticks = new(nticks,float)    ; Array to hold major tickmark values.
    imonth=0
    print("nticks="+nticks)
    do i=0,nticks-1
;      print("imonth="+imonth)
;      print("days_accum="+days_accum)
;      print("day_in_mon:"+day_in_mon(imonth))
      major_ticks(i)=days_accum
      days_accum     = days_accum + day_in_mon(imonth)
      imonth = imonth + 1
    end do

    res@tmXBMinorValues    = major_ticks(12:nticks-1:tk_per_year)
    res@tmXBValues         = major_ticks
    res@tmXBLabels         = new(nticks,string)
    res@tmXBLabels(:)     = ""

    print("nticks="+nticks)
;    res@tmXBLabels(7::tk_per_year)  = "Jan"
;    res@tmXBLabels(1::tk_per_year)  = "Y"+ispan(1,nyears+1,1)
 
    res@tmXBLabels(12::tk_per_year)  = ispan(1,nyears,1)
    res@tmXBLabels(6)  = "Years:"
    
  end if
  
  res@tmXBMode           = "Explicit"         ; Define own tick mark labels.

  res@tmXBMinorOutwardLengthF = 0.02    ; Increase length and size of
  res@tmXBMinorLengthF        = 0.01    ; minor tickmarks
  res@tmXBMinorThicknessF     = 2.0
  res@tmXBMajorThicknessF     = 2.0
  res@tmXBMajorOutwardLengthF = 0.01
  res@tmXBMajorLengthF        = 0.01
;===================================================================  
    
; Define our own contour levels.  
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnMinLevelValF=120
  res@cnMaxLevelValF=510
  res@cnLevelSpacingF     = 30.
  res@gsnRightString      = "O3 Column(DU)"
  res@lbLabelAngleF       = 60
  res@gsnLeftString       = "Control run"
;  plot = gsn_csm_lat_time(wks,O3_col_c(lat|:, time|:),res)      ; contour the variable
  res@gsnLeftString       = "Soot run"
;  plot = gsn_csm_lat_time(wks,O3_col_x(lat|:, time|:),res)      ; contour the variable
  
;  res@cnLevelSelectionMode = "ManualLevels" 
;  res@cnMinLevelValF=-330
;  res@cnMaxLevelValF=180
  res@cnLevelSelectionMode = "ExplicitLevels" 
  res@cnLevels          = (/-350,-300,-250,-200,-150,-100,-75, -50,-40,-30,-20,-10,-5,-2,0,10,20,40,80,120,160/) 
  res@cnLineThicknesses = (/  1 ,  1 ,  1 ,  1 ,  1 ,  1 , 1 ,  1 , 1 , 1 , 1 , 1 , 1, 1,2, 1, 1, 1, 1, 1 , 1/)
  res@cnMonoLineThickness = False
  res@cnLevelSpacingF     = 30.
  res@gsnLeftString       = "Soot run - Control run"
;  plot = gsn_csm_lat_time(wks,O3_col_diff(lat|:, time|:),res)      ; contour the variable

  res@cnLevelSelectionMode = "ExplicitLevels" 
  delete(res@cnLevels)
  res@cnLevels          = (/-60,-50,-40,-30,-20,-10,-5,-2,0,10,20,30,40,60,80,100,120/)
  delete(res@cnLineThicknesses)
  res@cnLineThicknesses = (/ 1 , 1 , 1 , 1 , 1 , 1 , 1, 1,4, 1, 1, 1, 1, 1, 1, 1 , 1 /)
  res@cnMonoLineThickness = False
  res@gsnLeftString       = "(Soot - Control)/Control"
  res@gsnRightString      = "O3 Column Change(%)"
  plot = gsn_csm_lat_time(wks,O3_col_dpct(lat|:, time|:),res)      ; contour the variable
end

undef("saveMeans")
procedure saveMeans(historyRun:string, controlRun:string)

begin
  
  saveMeanFiles(historyRun)
  
  saveMeanFiles(controlRun)

end

undef("plotO3Change")
procedure plotO3Change(historyRun:string, controlRun:string)
local T_x,T_c

begin
  print((/historyRun/))
;  saveMeanVars(historyRun,True)
;  saveMeanLOxFiles(historyRun)
  historyPattern = historyRun+"/"+historyRun+"_mean*.nc"
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  print((/historyNames/))
  O3_col_x=addfiles_GetVar(history, historyNames, "O3_col")
  O3_col_x=O3_col_x/2.687e16 ; Convert to Dobson Units
  O3_col_x@units="Dobson Units"
  print(dimsizes(O3_col_x))
 
;  T_x=File1->T
;  Soot=File1->SOOT01

  print((/controlRun/))
;  saveMeanVars(controlRun,False)
;  saveMeanLOxFiles(controlRun)
  controlPattern = controlRun+"/"+controlRun+"_mean*.nc"
  controlNames = systemfunc("ls " + controlPattern)
  control = addfiles(controlNames, "r")
  ListSetType(control, "cat")
  print((/controlNames/))
  O3_col_c=addfiles_GetVar(control, controlNames, "O3_col")
  O3_col_c=O3_col_c/2.687e16 ; Convert to Dobson Units
  O3_col_c@units="Dobson Units"
  print(dimsizes(O3_col_c))
  
  print("plotting ozone column...")
  plotO3Diff(O3_col_x, O3_col_c)
    
end
