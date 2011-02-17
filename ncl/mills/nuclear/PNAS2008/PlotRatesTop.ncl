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

undef("plotRatesTop")
procedure plotRatesTop(historyRun:string, controlRun:string)

begin
  print((/historyRun/))
  historyPattern = historyRun+"/"+historyRun+"_mean*.nc"
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  
  R_BrO_BrO_col_x = addfiles_GetVar(history, historyNames, "R_BrO_BrO_col")
  R_BrO_ClOa_col_x = addfiles_GetVar(history, historyNames, "R_BrO_ClOa_col")
  R_BrO_ClOb_col_x = addfiles_GetVar(history, historyNames, "R_BrO_ClOb_col")
  R_BrO_O_col_x = addfiles_GetVar(history, historyNames, "R_BrO_O_col")
  R_CH3O2_NO_col_x = addfiles_GetVar(history, historyNames, "R_CH3O2_NO_col")
  R_ClO_O_col_x = addfiles_GetVar(history, historyNames, "R_CLO_O_col")
  R_HO2_O_col_x = addfiles_GetVar(history, historyNames, "R_HO2_O_col")
  R_HO2_O3_col_x = addfiles_GetVar(history, historyNames, "R_HO2_O3_col")
  R_H_O3_col_x = addfiles_GetVar(history, historyNames, "R_H_O3_col")
  R_NO2_O_col_x = addfiles_GetVar(history, historyNames, "R_NO2_O_col")
  R_NO_HO2_col_x = addfiles_GetVar(history, historyNames, "R_NO_HO2_col")
  R_NO_N_col_x = addfiles_GetVar(history, historyNames, "R_NO_N_col")
  R_N_O2_col_x = addfiles_GetVar(history, historyNames, "R_N_O2_col")
  R_O1D_H2O_col_x = addfiles_GetVar(history, historyNames, "R_O1D_H2O_col")
  R_O2_O_col_x = addfiles_GetVar(history, historyNames, "R_O2_O_col")
  R_O3_O_col_x = addfiles_GetVar(history, historyNames, "R_O3_O_col")
  R_OH_O_col_x = addfiles_GetVar(history, historyNames, "R_OH_O_col")
  R_O_O_col_x = addfiles_GetVar(history, historyNames, "R_O_O_col")

  J_Cl2O2_col_x = addfiles_GetVar(history, historyNames, "J_Cl2O2_col")
  J_HOBr_col_x = addfiles_GetVar(history, historyNames, "J_HOBr_col")
  J_HOCl_col_x = addfiles_GetVar(history, historyNames, "J_HOCl_col")
  J_NO_col_x = addfiles_GetVar(history, historyNames, "J_NO_col")
  J_NO3_b_col_x = addfiles_GetVar(history, historyNames, "J_NO3_b_col")
  J_O2_a_col_x = addfiles_GetVar(history, historyNames, "J_O2_a_col")
  J_O2_b_col_x = addfiles_GetVar(history, historyNames, "J_O2_b_col")

  O3_col_x=addfiles_GetVar(history, historyNames, "O3_col")
  O3_col_x=O3_col_x/2.687e16 ; Convert to Dobson Units
  O3_col_x@units="Dobson Units"

  print((/controlRun/))
;  saveMeanVars(controlRun,False)
;  saveMeanLOxFiles(controlRun)
  controlPattern = controlRun+"/"+controlRun+"_mean*.nc"
  controlNames = systemfunc("ls " + controlPattern)
  control = addfiles(controlNames, "r")
  ListSetType(control, "cat")
;--------------------------------------------------------- 
;RE: Odd Oxygen Loss Rates (molecules cm-3 sec-1)
;O + O3
  R_O3_O_col_c = addfiles_GetVar(control, controlNames, "R_O3_O_col")
;NO2 + O
  R_NO2_O_col_c = addfiles_GetVar(control, controlNames, "R_NO2_O_col")
;JNO3_b [NO3]; NO + O2 branch
  J_NO3_b_col_c = addfiles_GetVar(control, controlNames, "J_NO3_b_col")
;HO2 + O
  R_HO2_O_col_c = addfiles_GetVar(control, controlNames, "R_HO2_O_col")
;HO2 + O3
  R_HO2_O3_col_c = addfiles_GetVar(control, controlNames, "R_HO2_O3_col")
;H + O3
  R_H_O3_col_c = addfiles_GetVar(control, controlNames, "R_H_O3_col")
;OH + O
  R_OH_O_col_c = addfiles_GetVar(control, controlNames, "R_OH_O_col")
;ClO + O
  R_ClO_O_col_c = addfiles_GetVar(control, controlNames, "R_CLO_O_col")
;ClO + BrO = Br + Cl + O2
  R_BrO_ClOa_col_c = addfiles_GetVar(control, controlNames, "R_BrO_ClOa_col")
;ClO + BrO = BrCl + O2
  R_BrO_ClOb_col_c = addfiles_GetVar(control, controlNames, "R_BrO_ClOb_col")
;JCl2O2 [Cl2O2]
  J_Cl2O2_col_c = addfiles_GetVar(control, controlNames, "J_Cl2O2_col")
;JHOCl [HOCl]
  J_HOCl_col_c = addfiles_GetVar(control, controlNames, "J_HOCl_col")
;BrO + O
  R_BrO_O_col_c = addfiles_GetVar(control, controlNames, "R_BrO_O_col")
;BrO + BrO
  R_BrO_BrO_col_c = addfiles_GetVar(control, controlNames, "R_BrO_BrO_col")
;JHOBr [HOBr]
  J_HOBr_col_c = addfiles_GetVar(control, controlNames, "J_HOBr_col")
;---------------------------------------------------------
;RE: Production (molecules cm-3 sec-1)
;JO2_a [O2] -> O + O(1D)
  J_O2_a_col_c = addfiles_GetVar(control, controlNames, "J_O2_a_col")
;JO2_b [O2] -> 2O
  J_O2_b_col_c = addfiles_GetVar(control, controlNames, "J_O2_b_col")
;HO2 + NO
  R_NO_HO2_col_c = addfiles_GetVar(control, controlNames, "R_NO_HO2_col")
;CH3O2 + NO
  R_CH3O2_NO_col_c = addfiles_GetVar(control, controlNames, "R_CH3O2_NO_col")
;---------------------------------------------------------
;RE: Other interesting rates (molecules cm-3 sec-1)
;O + O2 + M
  R_O2_O_col_c = addfiles_GetVar(control, controlNames, "R_O2_O_col")
;JO3_a[O3]
;JO3_b[O3]
;H2O + O1D
  R_O1D_H2O_col_c = addfiles_GetVar(control, controlNames, "R_O1D_H2O_col")
;N + O2
  R_N_O2_col_c = addfiles_GetVar(control, controlNames, "R_N_O2_col")
;NO + N
  R_NO_N_col_c = addfiles_GetVar(control, controlNames, "R_NO_N_col")
;JNO[NO]
  J_NO_col_c = addfiles_GetVar(control, controlNames, "J_NO_col")

  R_O_O_col_c = addfiles_GetVar(control, controlNames, "R_O_O_col")
  
  O3_col_c=addfiles_GetVar(control, controlNames, "O3_col")
  O3_col_c=O3_col_c/2.687e16 ; Convert to Dobson Units
  O3_col_c@units="Dobson Units"
    
  time=R_O3_O_col_x&time
  lat=R_O3_O_col_x&lat
 
  ; Create a panel plot of the loss rates.
  print("")
  j=-1
  res = True
  resP = True
  if (j.eq.-1) then
    print(" Plotting global mean...")
    wks = gsn_open_wks("pdf", "RatesTop_GlobAvg")
    resP@txString = "Global Mean"
  else
    print(" Plotting at latitude "+lat(j)+"...")
    wks = gsn_open_wks("pdf", "RatesTop_"+lat(j))
    resP@txString = "Latitude: "+lat(j)
  end if
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
    
  res@gsnMaximize           = True
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.4           ; change aspect ratio of plot
  res@vpHeightF           = 0.55

;  res@gsnDraw = False
;  res@gsnFrame = False
  
  res@trYMinF = 0
;  res@trYMaxF = 5.0
  
;  res@tmYLFormat = "0*+^e"
;  res@tmYLPrecision = 3
  res@pmLegendDisplayMode = "Always"
  res@pmLegendWidthF = 0.45
  res@lgLabelFontHeightF = 0.025
  res@pmLegendHeightF = 0.28
  
;------Page 1--------  

;  nPanel=3
;  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
;  plot = new(nPanel, graphic)
    
  res@gsnLeftString = "Control: Ox Column Production & Loss"
  res@gsnRightString = "molec/cm~S~2~N~/s"
  data= new((/6,dimsizes(time)/),float)
  data(0,:)=LOxLatMean(J_O2_b_col_c,j)
  data(1,:)=LOxLatMean(J_O2_b_col_x,j)
;2 x O + O3
  data(2,:)=2*LOxLatMean(R_O3_O_col_c,j)
  data(3,:)=2*LOxLatMean(R_O3_O_col_x,j)
;2 x NO2 + O
  data(4,:)=2*LOxLatMean(R_NO2_O_col_c,j)
  data(5,:)=2*LOxLatMean(R_NO2_O_col_x,j)
  
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

  res@xyLineThicknessF = 3.0
  res@xyLineColors      = (/"red","red","black","black","blue","blue"/); change line color
  res@xyDashPatterns  = (/1,0,1,0,1,0/)
  res@trXMaxF = time(dimsizes(time)-1)
  
  res@xyExplicitLegendLabels=(/"O~B~2~N~ + hv -> O(~S~1~N~D) + O control", "O~B~2~N~ + hv -> O(~S~1~N~D) + O soot", "2 x O~B~3~N~ + O control", "2 x O~B~3~N~ + O soot", "2 x NO~B~2~N~ + O control", "2 x NO~B~2~N~ + O soot" /)
  plot = gsn_csm_xy(wks, time, data, res)  
  
  delete(wks) 
    
end
