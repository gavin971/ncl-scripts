load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("LOxLatMean")
function LOxLatMean(Rate[*][*]:float,j:integer)

local pi,wgt
begin
  pi=acos(0.0)*2.0
  if (j.eq.-1) then
    Rate!1 = "lat"
    lat=Rate&lat
    wgt=dble2flt(cos(lat*pi/180.))
    LOx=dim_avg_wgt_Wrap(Rate,wgt,0)
    delete(wgt)
  else
    LOx=Rate(:,j)
  end if
  delete(pi)
return(LOx)
end


undef("plotO3Depletion")
procedure plotO3Depletion(historyRun:string, controlRun:string)

begin
  print((/historyRun/))
  historyPattern = historyRun+"/"+historyRun+".mean.h3.*.nc"
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
 
  controlPattern = controlRun+"/"+controlRun+".mean.h3.*.nc"
  controlNames = systemfunc("ls " + controlPattern)
  control = addfiles(controlNames, "r")
  ListSetType(control, "cat")

  O3_col_x=addfiles_GetVar(history, historyNames, "O3_col")
  O3_col_c=addfiles_GetVar(control, controlNames, "O3_col")
  print("O3_col_x:"+dimsizes(O3_col_x))
  print("O3_col_c:"+dimsizes(O3_col_c))
  dc=dimsizes(O3_col_c)
  dx=dimsizes(O3_col_x)
  d1=min((/dc(0),dx(0)/))
  
; smooth
  O3x=O3_col_x
  O3c=O3_col_c
  do i=5,d1-5
    O3_col_x(lat|:,time|i)=dim_avg(O3x(lat|:,time|i-5:i+4))
    O3_col_c(lat|:,time|i)=dim_avg(O3c(lat|:,time|i-5:i+4))
  end do
 
  O3_col_pct_loss= 100. * (O3_col_x(0:d1-1,:) - O3_col_c(0:d1-1,:)) / O3_col_c(0:d1-1,:)
  time=O3_col_c&time(0:d1-1)
  lat=O3_col_x&lat
  O3_col_pct_loss!1="lat"
  O3_col_pct_loss&lat=lat 
     
  wks = gsn_open_wks("pdf", "O3Depletion.h3.smooth")

  res = True
  resP = True
  resP@txString = "O3 Column % Change"
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
    
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.4           ; change aspect ratio of plot
  res@vpHeightF           = 0.55

;  res@gsnDraw = False
;  res@gsnFrame = False
  
;  res@trYMinF = 0
;  res@trYMaxF = 5.0
  
;  res@tmYLFormat = "0*+^e"
;  res@tmYLPrecision = 3
  res@pmLegendDisplayMode = "Always"
  nPanel=1
  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)
 
;===================================================================
; Do special labels on X axis. There will be 6 tickmarks per each
; month, with a month label appearing in the middle.  Here, we are
; treating minor tickmarks like major tickmarks, and vice versa, since
; only the major tickmarks can actually be labeled.
;
  year=1995
  yearlast=1997
  nyears=yearlast-year
  nmonths =12*nyears+2
  
  months  = new(nmonths,integer)
  years  = new(nmonths,integer)
  
  month1=1 ; Starts in January
  
;  months(0:11-month1) = ispan(month1,12,1)
;  years(0:11-month1) = year
  
;  months(0)=12
  
  do i=0,nmonths-3,12
    print(i)
    year=year+1
    years(i:i+11)=year
    months(i:i+11) =ispan(1,12,1) 
  end do
  months(nmonths-2:nmonths-1)=ispan(1,2,1)
  years(nmonths-2:nmonths-1)=yearlast
  
  day_in_mon  = new(nmonths,float)           ; Number of days in each month
  do i=0,nmonths-1
    day_in_mon(i) = days_in_month(years(i),months(i))*1.
  end do

  imonth      = 0
  days_accum  = 31.
  if (nmonths .lt. 60) then
    tk_per_mon  =4                                     ; Tickmarks per month.
    nticks       = nmonths * tk_per_mon
    major_ticks = new(nticks,float)    ; Array to hold major tickmark values.
    ii5         = ispan(1,tk_per_mon-1,1)*1.
    do i=2,nticks-3,tk_per_mon
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
    print("months: "+months)
    
   ; Label mid tickmark with the month label.
    res@tmXBLabels(0::tk_per_mon)  = monthnames(months-1)
;    print("tmXBLabels:"+res@tmXBLabels)
    
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

  res@tmYLMode           = "Manual" 
  res@tmYLTickStartF     = -70
  res@tmYLTickSpacingF        = 10
  res@tmYRMode           = "Manual"     
  res@tmYRTickStartF     = -70
  res@tmYRTickSpacingF   = 10
  res@tiYAxisString      = "Ozone column % change"
  
  res@pmLegendSide           = "Top"               ; Change location of 
;  res@pmLegendParallelPosF   = .45                  ; move units right
  res@pmLegendParallelPosF   = .5                  ; move units right
  res@pmLegendOrthogonalPosF = -0.9                ; move units down for 5Tg case
;  res@pmLegendOrthogonalPosF = -1.0                ; move units down for 1Tg case
;  res@lgLabelAngleF          = 290
;  res@lgOrientation = "horizontal"                 ; swap the orientation of the legend
  res@pmLegendWidthF         = 0.25                ; Change width and
  res@pmLegendHeightF        = 0.27                ; height of legend.
  res@lgPerimOn              = False               ; turn off box around
  res@lgLabelFontHeightF     = .035                 ; label font height
 
  res@tmYRLabelsOn           = True                ; Turn on right axis tick labels
  
  
  data= new((/5,dimsizes(time)/),float) 
  data(0,:)=O3_col_pct_loss(:,{-86}) ; 1Tg case
  data(1,:)=O3_col_pct_loss(:,{-46})
  data(2,:)=O3_col_pct_loss(:,{46})
  data(3,:)=O3_col_pct_loss(:,{86})
  data(4,:)=LOxLatMean(O3_col_pct_loss,-1)
  
; smooth
;  data1=data
;  do i=5,dimsizes(time)-5
;    data (:,i)=dim_avg(data1(:,i-5:i+4))
;  end do
  
  res@xyLineThicknesses = (/3,3,3,3,6/)
  res@xyMonoLineThickness = False
  res@xyLineColors      = (/"blue","magenta","DarkOliveGreen","red","black"/); change line color
  res@xyDashPatterns=(/0,0,0,0,0/)
  res@xyExplicitLegendLabels=(/"86S", "46S", "46N", "86N","global mean"/)
  res@trXMinF = 365.+15.
  res@trXMaxF = 365.*2.+15.
  res@trYMaxF = 15
  res@trYMinF = -45 
  res@xyLineThicknessF = 3.0
  res@tmXBLabelFontHeightF = 0.02
  plot(0) = gsn_csm_xy(wks, time, data, res)

;  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(plot)

end

;Main

plotO3Depletion("geo2", "nongeo2")
