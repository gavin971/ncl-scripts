load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "tests/getHistoryPath.ncl"

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

undef("plotRatesDiff")
procedure plotRatesDiff(historyRun:string, controlRun:string, j:integer)

begin
  print((/historyRun/))
  historyPath = getHistoryPath(historyRun)
  print((/historyPath/))
  historyPattern = historyPath+historyRun+".mean.h0.*.nc"
  historyNames = systemfunc("ls " + historyPattern)
  print((/historyNames/))
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
  controlPath = getHistoryPath(controlRun)
  print((/controlPath/))
;  saveMeanVars(controlRun,False)
;  saveMeanLOxFiles(controlRun)
  controlPattern = controlPath+controlRun+".mean.h0*.nc"
  controlNames = systemfunc("ls " + controlPattern)
  print((/controlNames/))
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
 
  print("")
  res = True
  if (j.eq.-1) then
    print(" Plotting global mean...")
    wks = gsn_open_wks("pdf", "RatesDiff.h0_GlobAvg")
    txString = "Global Mean"
  else
    print(" Plotting at latitude "+lat(j)+"...")
    wks = gsn_open_wks("pdf", "RatesDiff.h0_"+lat(j))
    txString = "Latitude: "+lat(j)
  end if
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
    
  res@gsnMaximize           = True
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.9

;  res@gsnDraw = False
;  res@gsnFrame = False
  
;  res@trYMinF = 0
;  res@trYMaxF = 5.0
  
;  res@tmYLFormat = "0*+^e"
;  res@tmYLPrecision = 3
;  res@pmLegendDisplayMode = "Always"
  res@pmLegendWidthF = 0.45
  res@lgLabelFontHeightF = 0.025
  res@pmLegendHeightF = 0.19
  res@pmLegendSide           = "Top"               ; Change location of
  res@pmLegendParallelPosF   = .45                 ; move units right
  res@pmLegendOrthogonalPosF = -1.05               ; move units down
  res@lgBoxLineThicknessF = 0.0                    ; remove box
  res@lgBoxLineColor="White"                       ; remove box
  
;------Page 1--------  

;  nPanel=3
;  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
;  plot = new(nPanel, graphic)
    
  res@gsnCenterString = txString+" Ox Column Loss Rates (Geoeng - Control)"
;  res@gsnRightString = "molec/cm~S~2~N~/s"
  dc=dimsizes(R_O3_O_col_c)
  dx=dimsizes(R_O3_O_col_x)
  d1=min((/dc(0),dx(0)/))
;  d0=365+31+15
  d0=0
;  d1=dx(0)
  data= new((/5,d1-d0/),float)

; ClOx
  data_x=LOxLatMean(R_ClO_O_col_x(d0:d1-1,:),j)+ LOxLatMean(J_Cl2O2_col_x(d0:d1-1,:),j)
  if (d1.gt.dc(0)) then
    data_c=data_x
    data_c(0:dc(d0)-1)=LOxLatMean(R_ClO_O_col_c(d0:dc(0)-1,:),j)+ \
                       LOxLatMean(J_Cl2O2_col_c(d0:dc(0)-1,:),j)
    data_c(dc(0):d1-1)=data_c(0:d1-dc(0)-1)
  else
    data_c=LOxLatMean(R_ClO_O_col_c(d0:d1-1,:),j)+ LOxLatMean(J_Cl2O2_col_c(d0:d1-1,:),j)
  end if
  data(0,:)=2e-8*(data_x-data_c)
  print("====================ClOx=============================")
  print((/data_x/)+" | "+(/data_c/)+" | "+(/data(0,:)/))
  print("=====================================================")

; Ox    
  data_x=LOxLatMean(R_O3_O_col_x(d0:d1-1,:),j)
  if (d1.gt.dc(0)) then
    data_c(d0:dc(0)-1)=LOxLatMean(R_O3_O_col_c(d0:dc(0)-1,:),j)
    data_c(dc(0):d1-1)=data_c(d0:d1-dc(0)-1)
  else
    data_c=LOxLatMean(R_O3_O_col_c(d0:d1-1,:),j)
  end if
  data(1,:)=2e-8*(data_x-data_c)
  print("=====================Ox==============================")
  print((/data_x/)+" | "+(/data_c/)+" | "+(/data(1,:)/))
  print("=====================================================")

; NOx   
  data_x=LOxLatMean(R_NO2_O_col_x(d0:d1-1,:),j)
  if (d1.gt.dc(0)) then
    data_c(d0:dc(0)-1)=LOxLatMean(R_NO2_O_col_c(d0:dc(0)-1,:),j)
    data_c(dc(0):d1-1)=data_c(d0:d1-dc(0)-1)
  else
    data_c=LOxLatMean(R_NO2_O_col_c(d0:d1-1,:),j)
  end if
  data(2,:)=2e-8*(data_x-data_c)
  print("====================NOx=============================")
  print((/data_x/)+" | "+(/data_c/)+" | "+(/data(2,:)/))
  print("=====================================================")

; HOx  
  data_x= LOxLatMean(R_OH_O_col_x(d0:d1-1,:),j)+ LOxLatMean(R_HO2_O_col_x(d0:d1-1,:),j)+ \
          LOxLatMean(R_HO2_O3_col_x(d0:d1-1,:),j)
  if (d1.gt.dc(0)) then
    data_c(d0:dc(0)-1)=LOxLatMean(R_OH_O_col_c(d0:dc(0)-1,:),j)+ \
                       LOxLatMean(R_HO2_O_col_c(d0:dc(0)-1,:),j)+ \
                       LOxLatMean(R_HO2_O3_col_c(d0:dc(0)-1,:),j)
    data_c(dc(0):d1-1)=data_c(d0:d1-dc(0)-1)
  else
    data_c= LOxLatMean(R_OH_O_col_c(d0:d1-1,:),j)+ LOxLatMean(R_HO2_O_col_c(d0:d1-1,:),j)+ \
            LOxLatMean(R_HO2_O3_col_c(d0:d1-1,:),j)
  end if
  data(3,:)=2e-8*(data_x-data_c)
  print("====================HOx=============================")
  print((/data_x/)+" | "+(/data_c/)+" | "+(/data(3,:)/))
  print("=====================================================")

; Brox             
  data_x=LOxLatMean(J_NO3_b_col_x(d0:d1-1,:),j)+ LOxLatMean(R_BrO_ClOa_col_x(d0:d1-1,:),j)+ \
         LOxLatMean(J_HOBr_col_x(d0:d1-1,:),j)+ LOxLatMean(R_BrO_ClOb_col_x(d0:d1-1,:),j)+ \
         LOxLatMean(R_BrO_O_col_x(d0:d1-1,:),j)+ LOxLatMean(R_BrO_BrO_col_x(d0:d1-1,:),j)
  if (d1.gt.dc(0)) then
    data_c(d0:dc(0)-1)=LOxLatMean(J_NO3_b_col_c(d0:dc(0)-1,:),j)+ LOxLatMean(R_BrO_ClOa_col_c(d0:dc(0)-1,:),j)+ \
           LOxLatMean(J_HOBr_col_c(d0:dc(0)-1,:),j)+ LOxLatMean(R_BrO_ClOb_col_c(d0:dc(0)-1,:),j)+ \
           LOxLatMean(R_BrO_O_col_c(d0:dc(0)-1,:),j)+ LOxLatMean(R_BrO_BrO_col_c(d0:dc(0)-1,:),j)
    data_c(dc(0):d1-1)=data_c(d0:d1-dc(0)-1)
  else
    data_c= LOxLatMean(J_NO3_b_col_c(d0:d1-1,:),j)+ LOxLatMean(R_BrO_ClOa_col_c(d0:d1-1,:),j)+ \
         LOxLatMean(J_HOBr_col_c(d0:d1-1,:),j)+ LOxLatMean(R_BrO_ClOb_col_c(d0:d1-1,:),j)+ \
         LOxLatMean(R_BrO_O_col_c(d0:d1-1,:),j)+ LOxLatMean(R_BrO_BrO_col_c(d0:d1-1,:),j)
  end if    
  data(4,:)=2e-8*(data_x-data_c)
  print("====================BrOx=============================")
  print((/data_x/)+" | "+(/data_c/)+" | "+(/data(4,:)/))
  print("=====================================================")
  
;===================================================================
; Do special labels on X axis. There will be 6 tickmarks per each
; month, with a month label appearing in the middle.  Here, we are
; treating minor tickmarks like major tickmarks, and vice versa, since
; only the major tickmarks can actually be labeled.
;
  nmonths =d1-d0+1
  year=1998
  nyears=floattoint(nmonths/12)+1  
  yearlast=year+nyears-1
  months  = new(nmonths,integer)
  years  = new(nmonths,integer)
  xmn = mod(nmonths,12)
  
;  month1=1 ; Starts in January
  
;  months(0:12-month1) = ispan(month1,12,1)
;  years(0:12-month1) = year
  
  do i=0,nmonths-xmn-1,12
    year=year+1
    years(i:i+11)=year
    months(i:i+11) =ispan(1,12,1) 
  end do

  if (xmn.gt.0) then
    months(nmonths-xmn:nmonths-1)=ispan(1,xmn,1)
    years(nmonths-xmn:nmonths-1)=yearlast
  end if
  
  day_in_mon  = new(nmonths,float)           ; Number of days in each month
  do i=0,nmonths-1
    day_in_mon(i) = days_in_month(years(i),months(i))*1.
  end do

  imonth      = 0
;  days_accum  = 31.
  days_accum  = dble2flt((/time(d0)/))
  print(days_accum)
  if (nmonths .lt. 60) then
    tk_per_mon  = 2                                ; Tickmarks per month.
    nticks       = nmonths * tk_per_mon
    major_ticks = new(nticks,float)    ; Array to hold major tickmark values.
    ii5         = ispan(1,tk_per_mon-1,1)*1.
    do i=0,nticks-1,tk_per_mon
      major_ticks(i) = days_accum
      print(i+" "+imonth+" "+days_accum+" "+nmonths)
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
;    print("labels:"+dimsizes(res@tmXBLabels(::tk_per_mon)))
    print("nticks/tk_per_mon="+nticks/tk_per_mon)

   ; Label mid tickmark with the month label.
    print((/months/))
    res@tmXBLabels(0::tk_per_mon)  = monthnames(months-1)
;    print("tmXBLabels:")
;    print((/res@tmXBValues/)+":"+(/res@tmXBLabels/))
    
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

  res@tmXBLabelFontHeightF        = 0.025
  res@xyLineThicknessF = 3.0
  res@xyLineColors      = (/"black","blue","red","green","brown"/); change line color
  res@xyDashPatterns  = (/0,0,0,0,0/)
;  res@xyLineColors      = (/"black","red","green"/); change line color
;  res@xyDashPatterns  = (/0,0,0/)
  res@trXMinF = time(d0)-15
  res@trXMaxF = time(d1-1)+15
  
;  res@xyExplicitLegendLabels=(/"O~B~2~N~ + hv -> O(~S~1~N~D) + O control", "O~B~2~N~ + hv -> O(~S~1~N~D) + O soot", "2 x O~B~3~N~ + O control", "2 x O~B~3~N~ + O soot", "2 x NO~B~2~N~ + O control", "2 x NO~B~2~N~ + O soot" /)
 ; res@xyExplicitLegendLabels=(/"2 x ClO + O", "2 x O~B~3~N~ + O", "2 x NO~B~2~N~ + O", "2 x Cl~B~2~N~O~B~2~N~ + hv" /)
  res@xyExplicitLegendLabels=(/"ClO~B~x~N~", "O~B~x~N~",  "NO~B~x~N~", "HO~B~x~N~", "BrO~B~x~N~"/)
;  res@xyExplicitLegendLabels=(/"ClO~B~x~N~", "NO~B~x~N~", "HO~B~x~N~"/)
  res@xyLineThicknessF=6
  res@tiYAxisString="Column loss (x10~S~8~N~ molec/cm~S~2~N~/s)"
  
  plot = gsn_csm_xy(wks, time(d0:d1-1), data, res)  
  
  delete(wks) 
    
end

; Main

;plotRatesDiff("geo6","nongeo5",2) ; 82S
;plotRatesDiff("geo6","nongeo5",-1) ; Global average
plotRatesDiff("geo6","nongeo5",43) ; 82N
