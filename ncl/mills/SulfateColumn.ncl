load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "tests/getHistoryPath.ncl"

undef("plotColumn")
procedure plotColumn(column[*][*]:float,filename:string) 
local wks,res,lev
begin 

  lat=column&lat
  
  xdim=dimsizes(column)
  print("xdim="+xdim)
  ntimes=xdim(0)
    
  wks = gsn_open_wks("pdf", filename)
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
  year=1995
  yearlast=1996
  nyears=yearlast-year
  nmonths =6+12*nyears
  
  months  = new(nmonths,integer)
  years  = new(nmonths,integer)
  months(0:11) = ispan(1,12,1)
  years(0:11) = year
  months(12:nmonths-1) = ispan(1,6,1)

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
;  res@cnLevelSelectionMode = "ManualLevels" 
;  res@cnMinLevelValF=120
;  res@cnMaxLevelValF=510
;  res@cnLevelSpacingF     = 30.
  res@gsnRightString      = "Sulfate Column(g/cm2)"
  res@lbLabelAngleF       = 60
;  res@gsnLeftString       = "Control run"
;  plot = gsn_csm_lat_time(wks,O3_col_c(lat|:, time|:),res)      ; contour the variable
;  res@gsnLeftString       = "Soot run"
;  plot = gsn_csm_lat_time(wks,O3_col_x(lat|:, time|:),res)      ; contour the variable
  
;  res@cnLevelSelectionMode = "ManualLevels" 
;  res@cnMinLevelValF=-330
;  res@cnMaxLevelValF=180
;  res@cnLevelSelectionMode = "ExplicitLevels" 
;  res@cnLevels          = (/-350,-300,-250,-200,-150,-100,-75, -50,-40,-30,-20,-10,-5,-2,0,10,20,40,80,120,160/) 
;  res@cnLineThicknesses = (/  1 ,  1 ,  1 ,  1 ,  1 ,  1 , 1 ,  1 , 1 , 1 , 1 , 1 , 1, 1,2, 1, 1, 1, 1, 1 , 1/)
;  res@cnMonoLineThickness = False
;  res@cnLevelSpacingF     = 30.
;  res@gsnLeftString       = "Soot run - Control run"
;  plot = gsn_csm_lat_time(wks,O3_col_diff(lat|:, time|:),res)      ; contour the variable

;  res@cnLevelSelectionMode = "ExplicitLevels" 
;  delete(res@cnLevels)
;  res@cnLevels          = (/-60,-50,-40,-30,-20,-10,-5,-2,0,10,20,30,40,60,80,100,120/)
;  delete(res@cnLineThicknesses)
;  res@cnLineThicknesses = (/ 1 , 1 , 1 , 1 , 1 , 1 , 1, 1,4, 1, 1, 1, 1, 1, 1, 1 , 1 /)
;  res@cnMonoLineThickness = False
;  res@gsnLeftString       = "(Soot - Control)/Control"
;  res@gsnRightString      = "O3 Column Change(%)"
  plot = gsn_csm_lat_time(wks,column(lat|:, time|:),res)      ; contour the variable
end

undef("plotSulfateColumn")
procedure plotSulfateColumn(historyRun:string)
local T_x,T_c

begin
  print((/historyRun/))
;  saveMeanVars(historyRun,True)
;  saveMeanLOxFiles(historyRun)
  historyPattern = getHistoryPath(historyRun)+historyRun+".mean.h3.nc"
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  print((/historyNames/))
  MIXSULMD_col=addfiles_GetVar(history, historyNames, "MIXSULMD_col")
  print(dimsizes(MIXSULMD_col))
  PURSULMD_col=addfiles_GetVar(history, historyNames, "PURSULMD_col")
  print(dimsizes(PURSULMD_col))
  TOTSULMD_col=MIXSULMD_col ; get metadata
  TOTSULMD_col=TOTSULMD_col+PURSULMD_col
   
  print("plotting sulfate column...")
  plotColumn(TOTSULMD_col,historyRun+".TOTSULMD.Column")
  plotColumn(MIXSULMD_col,historyRun+".MIXSULMD.Column")
  plotColumn(PURSULMD_col,historyRun+".PURSULMD.Column")
    
end

plotSulfateColumn("nogeoPurSpin1")
