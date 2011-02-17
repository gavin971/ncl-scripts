load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


undef("calcPress")
function calcPress(history:list, historyNames[*]:string)

local pi,pm,press,wgt
begin
  pi=acos(0.0)*2.0

  print("Reading surface pressure...")
  ps = addfiles_GetVar(history, historyNames, "PS")
  p0 = history[0]->P0
  hyam = history[0]->hyam
  hybm = history[0]->hybm
  time = addfiles_GetVar(history, historyNames, "time")
  lev = history[0]->lev
  lat = history[0]->lat
  lon = history[0]->lon
  
  print("Calculating pressure at midpoints...")
  pm = pres_hybrid_ccm(ps,p0,hyam,hybm)
  print(dimsizes(pm))
  pm!0 = "time"
  pm!1 = "lev"
  pm&lev = lev
  pm!2 = "lat"
  pm&lat = lat
  pm!3 = "lon"
  pm&lon = lon
  
  delete(p0)
  delete(ps)
  delete(hyam)
  delete(hybm)
  
  print("averaging pressure over time")
  press=dim_avg_Wrap(pm(lev|:, lat|:, lon|:, time|:))
  
  delete(pm)
  
  print("averaging pressure zonally")
  pm=dim_avg_Wrap(press(lev|:, lat|:, lon|:))
  
  delete(press)
  
  wgt=dble2flt(cos(lat*pi/180.))
  
  print("averaging pressure meridionally")
  pressure=dim_avg_wgt_Wrap(pm(lev|:, lat|:),wgt,0)
  
  delete(pm)
  
  pressure=pressure/100.
  
  pressure@units = "hPa"

  return(pressure)
end

undef("calcPdel")
function calcPdel(history:list, historyNames[*]:string)

local g,Mair,partPress,column,vmr,pdel
begin
 
  print("Reading surface pressure...")
  ps = addfiles_GetVar(history, historyNames, "PS")
  p0 = history[0]->P0
  hyai = history[0]->hyai
  hybi = history[0]->hybi
  time = addfiles_GetVar(history, historyNames, "time")
  lev = history[0]->lev
  lat = history[0]->lat
  lon = history[0]->lon
  
  print("Calculating pressure at interfaces...")
  pm = pres_hybrid_ccm(ps,p0,hyai,hybi)
    
  delete(p0)
  delete(ps)
  
  delete(hyai)
  delete(hybi)
  
  print("Calculating pressure change...")
  pdel=pm(:,1:dimsizes(lev),:,:)-pm(:,0:dimsizes(lev)-1,:,:)
  
  delete(pm)
  
  return(pdel)
    
end

undef("calcZdel")
function calcZdel(history:list, historyNames[*]:string)

local z3,zdel
begin
 
  print("Reading geopotential heights...")
  z3 = addfiles_GetVar(history, historyNames, "Z3")
  nlev=dimsizes(z3&lev)
  
  print("Calculating zdel...")
  zdel = z3
  zdel(:,0:nlev-2,:,:)=z3(:,0:nlev-2,:,:)-z3(:,1:nlev-1,:,:)
    
  delete(z3)
  delete(nlev)
  
  return(zdel)
    
end

undef("calcGlobMean")
function calcGlobMean(history:list, historyNames[*]:string, var:string, press[*]:float, wgt[*]:float)

local temp,tmean
begin

  print("reading "+var)
  temp=addfiles_GetVar(history, historyNames, var)
  print(dimsizes(temp))
  
  print("averaging zonally")
  tmean=dim_avg_Wrap(temp(time|:, lev|:, lat|:, lon|:))
  
  delete(temp)
  
  print("averaging meridionally")
  temp=dim_avg_wgt_Wrap(tmean(time|:, lev|:, lat|:),wgt,0)
  temp!0 = "time"
  temp!1 = "press"
  temp&press = press
  
  return(temp)

end


undef("plotTDiff")
procedure plotTDiff(T_x[*][*]:float, T_c[*][*]:float, Soot[*][*]:float, SootMass:float)
local T_diff,T_x,T_c,wks,res,lev
begin
  press=T_x&press
  print("pressure range:"+press(0)+" to "+press(dimsizes(press)-1))
  T_diff=T_x-T_c
  T_diff!0 = "time"
  T_diff!1 = "press"
  T_diff&press = press
  
  time=T_x&time

  wks = gsn_open_wks("pdf", "SootTemp")
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  res                     = True          ; plot mods desired
  resP = True

  res@gsnDraw = False
  res@gsnFrame = False
  
  nPanel=2
  res@vpHeightF             = 0.8 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)

  res@cnFillOn            = True          ; turn on color
;  res@lbLabelBarOn        = False         ; turn off color labelbar
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@gsnMaximize         = True          ; expand plot size to maximum
  resP@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  resP@gsnPaperOrientation = "landscape"
  resP@gsnPaperHeight = 12
  res@gsnPaperHeight = 12
  res@gsnPaperMargin = 0.0
  res@vpWidthF            = 1.55           ; change aspect ratio of plot
;  resP@vpHeightF            = 1.0           ; change aspect ratio of plot
  res@cnLevelSpacingF     = 10.
  res@gsnRightString      = "(K)"
  res@gsnLeftStringFontHeightF = 48
  res@gsnRightStringFontHeightF = 48
  res@tiYAxisString       = "Pressure (hPa)"
;  res@tiXAxisString       = "Days after May 15"
  res@cnLineLabelPlacementMode= "Computed"
  res@cnLineLabelDensityF = 5 
  
;===================================================================
; Do special labels on X axis. There will be 6 tickmarks per each
; month, with a month label appearing in the middle.  Here, we are
; treating minor tickmarks like major tickmarks, and vice versa, since
; only the major tickmarks can actually be labeled.
;
  nmonths =8+12+12+10
  months  = new(nmonths,integer)
  months(0:7) = (/ 5, 6, 7, 8, 9,10,11,12/)
  months(8:19) =(/ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12/) 
  months(20:31) =(/ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12/) 
  months(32:nmonths-1) =(/ 1, 2, 3, 4, 5, 6, 7, 8, 9,10/) 
  years  = new(nmonths,integer)
  years(0:7) = 2025
  years(8:19) =2026 
  years(20:31) =2027
  years(32:nmonths-1) =2028
  tk_per_mon  =4                                     ; Tickmarks per month.
;  nmonths     = dimsizes(months)
  ndays       = nmonths * tk_per_mon

  day_in_mon  = new(nmonths,float)           ; Number of days in each month
  do i=0,nmonths-1
    day_in_mon(i) = days_in_month(years(i),months(i))*1.
  end do

  major_ticks = new(ndays,float)    ; Array to hold major tickmark values.
  imonth      = 0
  days_accum  = 0.
  ii5         = ispan(1,tk_per_mon-1,1)*1.
  do i=0,ndays-1,tk_per_mon
    major_ticks(i) = days_accum
    major_ticks(i+1:i+tk_per_mon-1) = major_ticks(i) + ii5*(day_in_mon(imonth)/tk_per_mon)
    days_accum     = days_accum + day_in_mon(imonth)
    imonth = imonth + 1
  end do

  res@tmXBMode           = "Explicit"         ; Define own tick mark labels.
  res@tmXBMinorValues    = major_ticks(tk_per_mon/2:ndays-1:tk_per_mon)
  res@tmXBValues         = major_ticks
  res@tmXBLabels         = new(ndays,string)
  res@tmXBLabels(:)     = ""

  print("ndays="+ndays)
  monthnames=(/"J","F","M","A","M","J","J","A","S","O","N","D"/)
  print("monthnames:"+dimsizes(monthnames))
  print("labels:"+dimsizes(res@tmXBLabels(::tk_per_mon)))
  print("ndays/tk_per_mon="+ndays/tk_per_mon)
  
 ; Label mid tickmark with the month label.
  res@tmXBLabels(::tk_per_mon)  = monthnames(months-1)

  res@tmXBMinorOutwardLengthF = 0.02    ; Increase length and size of
  res@tmXBMinorLengthF        = 0.01    ; minor tickmarks
  res@tmXBMinorThicknessF     = 2.0
  res@tmXBMajorThicknessF     = 2.0
  res@tmXBMajorOutwardLengthF = 0.01
  res@tmXBMajorLengthF        = 0.01
;===================================================================  
 res@gsnLeftString       = "Total soot mass"
 res@gsnRightString       = "Tg"
 res@tiYAxisString       = "Soot"
 res@trYMinF = 0
 res@trYMaxF = 5.0
 res@xyLineThicknessF = 2.0
 res@tmYRLabelsOn           = True                ; Turn on right axis tick labels
 plot(1) = gsn_csm_xy(wks,time,SootMass/1e9,res)

  res@tmYLMode   = "Explicit" 
;  res@tmYLValues  = (/    1E-5 ,    3E-5   ,   1E-4  ,    3E-4   , 0.001 , 0.003 , 0.01 ,  0.03 , 0.1 , 0.3 , 1., 3., 10., 30., 100., 300./)
;  res@tmYLLabels  = (/"10~S~-5","3x10~S~-5","10~S~-4","3x10~S~-4","0.001","0.003","0.01", "0.03","0.1","0.3","1","3","10","30","100","300"/)
  res@tmYLValues  = (/    1E-5 ,   1E-4  , 0.001 , 0.01 , 0.1 , 1., 10., 100./)
  res@tmYLLabels  = (/"10~S~-5","10~S~-4","0.001","0.01","0.1","1","10","100"/)
;  res@tmLabelAutoStride = True  
 
;  res@gsnLeftString       = "Soot run, Temperature"
;  plot = gsn_csm_pres_hgt (wks,T_x(press|:, time|:),res)   ; contour the variable
;  res@gsnLeftString       = "Control run, Temperature"
;  plot = gsn_csm_pres_hgt(wks,T_c(press|:, time|:),res)      ; contour the variable
  
  res@gsnLeftString       = "Temperature, Soot - Control, global average"
  res@cnLineLabelInterval = 1 ; label every contour
;  res@cnMonoLevelFlag = True
;  res@cnLevelFlag = "LineAndLabel"
  res@cnLineLabelPlacementMode = "Computed"
;  res@cnLineDashSegLenF = 0.05
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnMinLevelValF=-20
  res@cnMaxLevelValF=100
  res@cnLevelSpacingF     = 10.
;  plot(1) = gsn_csm_pres_hgt(wks,T_diff(press|:, time|:),res)      ; contour the variable
  
; Define our own contour levels.  
  res@cnLevelSelectionMode = "ExplicitLevels" 
  res@cnLevels  = (/1,2,3,5,7,10,20,30,50,70,100/) 
;  res@cnLineLabelStrings = (/"1","2","3","5","7","10","20","30","50","70","100"/) 
  Soot=Soot*1e9
  res@gsnLeftString       = "Soot mmr, global average"
  res@gsnRightString      = "(10~S~-9~N~ kg/kg air)"
  res@tiYAxisString       = "Pressure (hPa)"
  delete(res@trYMinF)
  delete(res@trYMaxF)
  plot(0) = gsn_csm_pres_hgt(wks,Soot(press|:, time|:),res)      ; contour the variable

  
  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(plot)  
  
end


undef("plotSootTemp")
procedure plotSootTemp(historyRun:string, controlRun:string)
local T_x,T_c

begin
  print((/historyRun/))
  historyPattern = historyRun+"/"+historyRun+"_mean*.nc"
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  T_x=addfiles_GetVar(history, historyNames, "T")
  Soot=addfiles_GetVar(history, historyNames, "SOOT01")

  print((/controlRun/))
  controlPattern = controlRun+"/"+controlRun+"_mean*.nc"
  controlNames = systemfunc("ls " + controlPattern)
  control = addfiles(controlNames, "r")
  ListSetType(control, "cat")
  T_c=addfiles_GetVar(control, controlNames, "T")

  historyPattern2 = historyRun+"/SootMass*.nc"
  historyNames2 = systemfunc("ls " + historyPattern2)
  print(historyNames2)
  history2 = addfiles(historyNames2, "r")
  ListSetType(history, "cat")
  SootMass=addfiles_GetVar(history2, historyNames2, "SootMass")
   
  print("plotting Temperature...")
  plotTDiff(T_x,T_c,Soot,SootMass)
    
end
