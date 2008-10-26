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

undef("plotLOx")
procedure plotLOx(historyRun:string, controlRun:string, j:integer)

begin
  print((/historyRun/))
  historyPattern = historyRun+"/"+historyRun+"_mean199*.nc"
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
  controlPattern = controlRun+"/"+controlRun+"_mean199*.nc"
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
;  j=-1
  res = True
  resP = True
  if (j.eq.-1) then
    print(" Plotting global mean...")
    wks = gsn_open_wks("pdf", "LOx_GlobAvg")
    resP@txString = "Global Mean"
  else
    print(" Plotting at latitude "+lat(j)+"...")
    wks = gsn_open_wks("pdf", "LOx_"+lat(j))
    resP@txString = "Latitude: "+lat(j)
  end if
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  
  resP@vpWidthF              = 1.0
  resP@vpHeightF             = 1.0
  resP@gsnMaximize           = True
  
;  res@gsnMaximize           = True
  res@vpWidthF              = 1.0
  res@gsnDraw = False
  res@gsnFrame = False
  
  res@trYMinF = 0
;  res@trYMaxF = 5.0
  
;  res@tmYLFormat = "0*+^e"
;  res@tmYLPrecision = 3
  res@pmLegendDisplayMode = "Always"
  res@pmLegendWidthF = 0.4
;------Page 1--------  

  nPanel=3
  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)
  
  data= new((/2,dimsizes(time)/),float)
  
  res@gsnRightString = "DU"
  res@gsnLeftString = "O3 column"
  data= new((/2,dimsizes(time)/),float)
  print("time:")
  print(dimsizes(time))
  print("O3_col_x:")
  print(dimsizes(O3_col_x))
  print("O3_col_c:")
  print(dimsizes(O3_col_c))
  print("data:")
  print(dimsizes(data))
  
  data(0,:)=LOxLatMean(O3_col_x,j)
  data(1,:)=LOxLatMean(O3_col_c,j)
  res@xyExplicitLegendLabels=(/"Geoeng","Control"/)
  plot(0) = gsn_csm_xy(wks, time, data, res)  
  
  delete(data)
  res@gsnLeftString = "Control: Ox Column Production"
  res@gsnRightString = "molec/cm~S~2~N~/s"
  data= new((/4,dimsizes(time)/),float)
;  Rate=J_O2_a_col_c+J_O2_b_col_c
;  Rate!1="lat"
;  Rate&lat=lat
  data(0,:)=LOxLatMean(J_O2_a_col_c,j)
  data(1,:)=LOxLatMean(J_O2_b_col_c,j)
  data(2,:)=LOxLatMean(R_NO_HO2_col_c,j)
  data(3,:)=2*LOxLatMean(R_CH3O2_NO_col_c,j)
  res@xyLineThicknessF = 2.0
  res@xyLineColors      = (/"blue","red","green","black"/); change line color
  delete(res@xyExplicitLegendLabels)
  res@xyExplicitLegendLabels=(/"O~B~2~N~ + hv -> O(1D) + O", "O~B~2~N~ + hv -> 2O", "NO + HO~B~2~N~", "2 x CH3O2 + NO"/)
  plot(1) = gsn_csm_xy(wks, time, data, res)  
  
  delete(data)
  res@gsnLeftString = "Geoeng: Ox Column Production"
  res@gsnRightString = "molec/cm~S~2~N~/s"
  data= new((/4,dimsizes(time)/),float)
  data(0,:)=LOxLatMean(J_O2_a_col_x,j)
  data(1,:)=LOxLatMean(J_O2_b_col_x,j)
  data(2,:)=LOxLatMean(R_NO_HO2_col_x,j)
  data(3,:)=2*LOxLatMean(R_CH3O2_NO_col_x,j)
  res@xyLineThicknessF = 2.0
  res@pmLegendDisplayMode = "NoCreate"
  res@xyLineColors      = (/"blue","red","green","black"/); change line color
;  delete(res@xyExplicitLegendLabels)
;  res@xyExplicitLegendLabels=(/"O~B~2~N~+hv (a)", "O~B~2~N~+hv (b)", "NO+HO~B~2~N~", "CH~D~3~N~O~D~2~N~+NO"/)
  plot(2) = gsn_csm_xy(wks, time, data, res)  
  
  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(plot)

;------Page 2--------  
  res@pmLegendDisplayMode = "Always"
  
  nPanel=1
  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)
  
  delete(data)
  res@gsnLeftString = "Ox Column Loss"
  res@gsnRightString = "Geoeng/Control"
  data= new((/15,dimsizes(time)/),float)
;2 x O + O3
  data(0,:)=LOxLatMean(R_O3_O_col_x,j)/LOxLatMean(R_O3_O_col_c,j)
;2 x NO2 + O
  data(1,:)=LOxLatMean(R_NO2_O_col_x,j)/LOxLatMean(R_NO2_O_col_c,j)
;HO2 + O
  data(2,:)=LOxLatMean(R_HO2_O_col_x,j)/LOxLatMean(R_HO2_O_col_c,j)
;OH + O
  data(3,:)=LOxLatMean(R_OH_O_col_x,j)/LOxLatMean(R_OH_O_col_c,j)
;H + O3
  data(4,:)=LOxLatMean(R_H_O3_col_x,j)/LOxLatMean(R_H_O3_col_c,j)
;HO2 + O3
  data(5,:)=LOxLatMean(R_HO2_O3_col_x,j)/LOxLatMean(R_HO2_O3_col_c,j)
;2 x ClO + O
  data(6,:)=LOxLatMean(R_ClO_O_col_x,j)/LOxLatMean(R_ClO_O_col_c,j)
;HOCl + hv
  data(7,:)=LOxLatMean(J_HOCl_col_x,j)/(LOxLatMean(J_HOCl_col_c,j)+1e-30)
;BrO + O
  data(8,:)=LOxLatMean(R_BrO_O_col_x,j)/LOxLatMean(R_BrO_O_col_c,j)
;BrO + BrO
  data(9,:)=LOxLatMean(R_BrO_BrO_col_x,j)/LOxLatMean(R_BrO_BrO_col_c,j)
;ClO + BrO = BrCl
  data(10,:)=LOxLatMean(R_BrO_ClOa_col_x,j)/LOxLatMean(R_BrO_ClOa_col_c,j)
;ClO + BrO = Br + Cl
  data(11,:)=LOxLatMean(R_BrO_ClOb_col_x,j)/LOxLatMean(R_BrO_ClOb_col_c,j)
;JCl2O2 [Cl2O2]
  data(12,:)=LOxLatMean(J_Cl2O2_col_x,j)/(LOxLatMean(J_Cl2O2_col_c,j)+1e-30)
;JHOBr [HOBr]
  data(13,:)=LOxLatMean(J_HOBr_col_x,j)/(LOxLatMean(J_HOBr_col_c,j)+1e-30)
;NO3 + hv
  data(14,:)=LOxLatMean(J_NO3_b_col_x,j)/(LOxLatMean(J_NO3_b_col_c,j)+1e-30)
  
  
  res@xyLineThicknessF = 2.0
  delete(res@xyLineColors)
  res@xyLineColors      = (/"blue", "red", "green", "black", "yellow", "orange", "purple", "blue", "red", "green", "black", "yellow", "orange", "purple", "blue"/); change line color
  res@xyDashPatterns=(/0,0,0,0,0,0,0,1,1,1,1,1,1,1,2/)
  
  res@pmLegendDisplayMode = "Always"
  delete(res@xyExplicitLegendLabels)
  res@xyExplicitLegendLabels=(/"2 x O~B~3~N~ + O", "2 x NO~B~2~N~ + O", "HO2 + O", "OH + O", "H + O3", "HO2 + O3", "2 x ClO + O", "HOCl + hv", "2 x BrO + O", "2 x BrO + BrO", "2 x ClO + BrO -> Br + Cl + O2", "2 x ClO + BrO -> BrCl + O2", "2 x  Cl2O2+hv", "2 x HOBr+hv", "2 x NO3 + hv"/)
  plot(0) = gsn_csm_xy(wks, time, data, res) 
   
  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(plot)
  
;------Page 3--------  
  res@pmLegendDisplayMode = "Always"
  
  nPanel=1
  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)
  
  delete(data)
  res@gsnLeftString = "Ox Column Loss"
  res@gsnRightString = "Geoeng/Control"
  data= new((/14,dimsizes(time)/),float)
;2 x O + O3
  data(0,:)=LOxLatMean(R_O3_O_col_x,j)/(LOxLatMean(R_O3_O_col_c,j)+1e-30)
;2 x NO2 + O
  data(1,:)=LOxLatMean(R_NO2_O_col_x,j)/(LOxLatMean(R_NO2_O_col_c,j)+1e-30)
;HO2 + O
  data(2,:)=LOxLatMean(R_HO2_O_col_x,j)/(LOxLatMean(R_HO2_O_col_c,j)+1e-30)
;OH + O
  data(3,:)=LOxLatMean(R_OH_O_col_x,j)/(LOxLatMean(R_OH_O_col_c,j)+1e-30)
;H + O3
  data(4,:)=LOxLatMean(R_H_O3_col_x,j)/(LOxLatMean(R_H_O3_col_c,j)+1e-30)
;HO2 + O3
  data(5,:)=LOxLatMean(R_HO2_O3_col_x,j)/(LOxLatMean(R_HO2_O3_col_c,j)+1e-30)
;2 x ClO + O
  data(6,:)=LOxLatMean(R_ClO_O_col_x,j)/(LOxLatMean(R_ClO_O_col_c,j)+1e-30)
;HOCl + hv
  data(7,:)=LOxLatMean(J_HOCl_col_x,j)/(LOxLatMean(J_HOCl_col_c,j)+1e-30)
;BrO + O
  data(8,:)=LOxLatMean(R_BrO_O_col_x,j)/(LOxLatMean(R_BrO_O_col_c,j)+1e-30)
;BrO + BrO
  data(9,:)=LOxLatMean(R_BrO_BrO_col_x,j)/(LOxLatMean(R_BrO_BrO_col_c,j)+1e-30)
;ClO + BrO = BrCl
  data(10,:)=LOxLatMean(R_BrO_ClOa_col_x,j)/(LOxLatMean(R_BrO_ClOa_col_c,j)+1e-30)
;ClO + BrO = Br + Cl
  data(11,:)=LOxLatMean(R_BrO_ClOb_col_x,j)/(LOxLatMean(R_BrO_ClOb_col_c,j)+1e-30)
;JCl2O2 [Cl2O2]
  data(12,:)=LOxLatMean(J_Cl2O2_col_x,j)/(LOxLatMean(J_Cl2O2_col_c,j)+1e-30)
;JHOBr [HOBr]
  data(13,:)=LOxLatMean(J_HOBr_col_x,j)/(LOxLatMean(J_HOBr_col_c,j)+1e-30)
  
  
  res@xyLineThicknessF = 2.0
  delete(res@xyLineColors)
  res@xyLineColors      = (/"blue", "red", "green", "black", "yellow", "orange", "purple", "blue", "red", "green", "black", "yellow", "orange", "purple"/); change line color
  delete(res@xyDashPatterns)
  res@xyDashPatterns=(/0,0,0,0,0,0,0,1,1,1,1,1,1,1/)
  res@pmLegendDisplayMode = "Always"
  delete(res@xyExplicitLegendLabels)
  res@xyExplicitLegendLabels=(/"2 x O~B~3~N~ + O", "2 x NO~B~2~N~ + O", "HO2 + O", "OH + O", "H + O3", "HO2 + O3", "2 x ClO + O", "HOCl + hv", "2 x BrO + O", "2 x BrO + BrO", "2 x ClO + BrO -> Br + Cl + O2", "2 x ClO + BrO -> BrCl + O2", "2 x  Cl2O2+hv", "2 x HOBr+hv"/)
  res@trYMaxF = 4.0
  plot(0) = gsn_csm_xy(wks, time, data, res) 
   
  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(plot)
  delete(res@trYMaxF)
  
;------Page 4--------  
  res@pmLegendDisplayMode = "Always"
  
  nPanel=4
  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)
  
  delete(data)
  res@gsnLeftString = "Control: Ox Column Loss"
  res@gsnRightString = "molec/cm~S~2~N~/s"
  data= new((/4,dimsizes(time)/),float)
;2 x O + O3
  data(0,:)=2*LOxLatMean(R_O3_O_col_c,j)
;2 x NO2 + O
  data(1,:)=2*LOxLatMean(R_NO2_O_col_c,j)
;HO2 + O
  data(2,:)=LOxLatMean(R_HO2_O_col_c,j)
;OH + O
  data(3,:)=LOxLatMean(R_OH_O_col_c,j)
  res@xyLineThicknessF = 2.0
  delete(res@xyLineColors)
  res@xyLineColors      = (/"blue", "red", "green", "black"/); change line color
  res@pmLegendDisplayMode = "Always"
  delete(res@xyExplicitLegendLabels)
  res@xyExplicitLegendLabels=(/"2 x O~B~3~N~ + O", "2 x NO~B~2~N~ + O", "HO2 + O", "OH + O", "H + O3"/)
  plot(0) = gsn_csm_xy(wks, time, data, res)  
    
  delete(data)
  res@gsnLeftString = "Geoeng: Ox Column Loss"
  res@gsnRightString = "molec/cm~S~2~N~/s"
  data= new((/4,dimsizes(time)/),float)
;2 x O + O3
  data(0,:)=2*LOxLatMean(R_O3_O_col_x,j)
;2 x NO2 + O
  data(1,:)=2*LOxLatMean(R_NO2_O_col_x,j)
;HO2 + O
  data(2,:)=LOxLatMean(R_HO2_O_col_x,j)
;OH + O
  data(3,:)=LOxLatMean(R_OH_O_col_x,j)
;  res@xyLineThicknessF = 2.0
;  res@xyLineColors      = (/"blue", "red", "green", "black"/); change line color
  res@pmLegendDisplayMode = "NoCreate"
;  delete(res@xyExplicitLegendLabels)
;  res@xyExplicitLegendLabels=(/"2 x O~B~3~N~ + O", "2 x ClO + O", "2 x NO~B~2~N~ + O", "NO~D~3~N~ + hv", "HO2 + O", "HO2 + O3", "H + O3", "OH + O", "ClO + BrO -> BrCl", "ClO + BrO -> Br + Cl", "Cl2O2+hv", "BrO + O", "BrO + BrO", "HOBr+hv"/)
  plot(1) = gsn_csm_xy(wks, time, data, res)  
  
  delete(data)
  res@gsnLeftString = "Control: Ox Column Loss"
  res@gsnRightString = "molec/cm~S~2~N~/s"
  data= new((/4,dimsizes(time)/),float)
;H + O3
  data(0,:)=LOxLatMean(R_H_O3_col_c,j)
;HO2 + O3
  data(1,:)=LOxLatMean(R_HO2_O3_col_c,j)
;2 x ClO + O
  data(2,:)=2*LOxLatMean(R_ClO_O_col_c,j)
;NO3 + hv
  data(3,:)=2*LOxLatMean(J_NO3_b_col_c,j)
  res@xyLineThicknessF = 2.0
  delete(res@xyLineColors)
  res@xyLineColors      = (/"blue", "red", "green", "black", "orange", "DarkSlateGray", "magenta", "DarkOliveGreen"/); change line color
  res@pmLegendDisplayMode = "Always"
  delete(res@xyExplicitLegendLabels)
  res@xyExplicitLegendLabels=(/"H + O3", "HO2 + O3", "2 x ClO + O", "2 x NO3 + hv"/)
  plot(2) = gsn_csm_xy(wks, time, data, res)  
    
  delete(data)
  res@gsnLeftString = "Geoeng: Ox Column Loss"
  res@gsnRightString = "molec/cm~S~2~N~/s"
  data= new((/4,dimsizes(time)/),float)
;H + O3
  data(0,:)=LOxLatMean(R_H_O3_col_x,j)
;HO2 + O3
  data(1,:)=LOxLatMean(R_HO2_O3_col_x,j)
;2 x ClO + O
  data(2,:)=2*LOxLatMean(R_ClO_O_col_x,j)
;NO3 + hv
  data(3,:)=2*LOxLatMean(J_NO3_b_col_x,j)
;  res@xyLineThicknessF = 2.0
;  res@xyLineColors      = (/"blue", "red", "green", "black", "orange", "DarkSlateGray", "magenta", "DarkOliveGreen"/); change line color
  res@pmLegendDisplayMode = "NoCreate"
;  delete(res@xyExplicitLegendLabels)
;  res@xyExplicitLegendLabels=(/"2 x O~B~3~N~ + O", "2 x ClO + O", "2 x NO~B~2~N~ + O", "NO~D~3~N~ + hv", "HO2 + O", "HO2 + O3", "H + O3", "OH + O", "ClO + BrO -> BrCl", "ClO + BrO -> Br + Cl", "Cl2O2+hv", "BrO + O", "BrO + BrO", "HOBr+hv"/)
  plot(3) = gsn_csm_xy(wks, time, data, res)  
  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(plot)

;------Page 5--------  
  res@pmLegendDisplayMode = "Always"
  
  nPanel=4
  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)

  delete(data)
  res@gsnLeftString = "Control: Ox Column Loss"
  res@gsnRightString = "molec/cm~S~2~N~/s"
  data= new((/3,dimsizes(time)/),float)
;HOCl + hv
  data(0,:)=LOxLatMean(J_HOCl_col_c,j)
;BrO + O
  data(1,:)=2*LOxLatMean(R_BrO_O_col_c,j)
;BrO + BrO
  data(2,:)=2*LOxLatMean(R_BrO_BrO_col_c,j)
  res@xyLineThicknessF = 2.0
  delete(res@xyLineColors)
  res@xyLineColors      = (/"blue", "red", "green", "black", "orange", "DarkSlateGray", "magenta"/); change line color
  res@pmLegendDisplayMode = "Always"
  delete(res@xyExplicitLegendLabels)
  res@xyExplicitLegendLabels=(/"HOCl + hv", "2 x BrO + O", "2 x BrO + BrO"/)
  plot(0) = gsn_csm_xy(wks, time, data, res)  
    
  delete(data)
  res@gsnLeftString = "Geoeng: Ox Column Loss"
  res@gsnRightString = "molec/cm~S~2~N~/s"
  data= new((/3,dimsizes(time)/),float)
;HOCl + hv
  data(0,:)=LOxLatMean(J_HOCl_col_x,j)
;BrO + O
  data(1,:)=2*LOxLatMean(R_BrO_O_col_x,j)
;BrO + BrO
  data(2,:)=2*LOxLatMean(R_BrO_BrO_col_x,j)
  res@xyLineThicknessF = 2.0
;  res@xyLineColors      = (/"blue", "red", "green", "black", "orange", "DarkSlateGray", "magenta"/); change line color
  res@pmLegendDisplayMode = "NoCreate"
;  delete(res@xyExplicitLegendLabels)
;  res@xyExplicitLegendLabels=(/"2 x O~B~3~N~ + O", "2 x ClO + O", "2 x NO~B~2~N~ + O", "NO~D~3~N~ + hv", "HO2 + O", "HO2 + O3", "H + O3", "OH + O", "ClO + BrO -> BrCl", "ClO + BrO -> Br + Cl", "Cl2O2+hv", "BrO + O", "BrO + BrO", "HOBr+hv"/)
  plot(1) = gsn_csm_xy(wks, time, data, res)  
  
  delete(data)
  res@gsnLeftString = "Control: Ox Column Loss"
  res@gsnRightString = "molec/cm~S~2~N~/s"
  data= new((/4,dimsizes(time)/),float)
;ClO + BrO = Br + Cl + O2
  data(0,:)=2*LOxLatMean(R_BrO_ClOa_col_c,j)
;ClO + BrO = BrCl + O2
  data(1,:)=2*LOxLatMean(R_BrO_ClOb_col_c,j)
;JCl2O2 [Cl2O2]
  data(2,:)=2*LOxLatMean(J_Cl2O2_col_c,j)
;JHOBr [HOBr]
  data(3,:)=2*LOxLatMean(J_HOBr_col_c,j)
  res@xyLineThicknessF = 2.0
  delete(res@xyLineColors)
  res@xyLineColors      = (/"blue", "red", "green", "black"/); change line color
  res@pmLegendDisplayMode = "Always"
  delete(res@xyExplicitLegendLabels)
  res@xyExplicitLegendLabels=(/"2 x ClO + BrO -> Br + Cl + O2", "2 x ClO + BrO -> BrCl + O2", "2 x  Cl2O2+hv", "2 x HOBr+hv"/)
  plot(2) = gsn_csm_xy(wks, time, data, res)  
    
  delete(data)
  res@gsnLeftString = "Geoeng: Ox Column Loss"
  res@gsnRightString = "molec/cm~S~2~N~/s"
  data= new((/4,dimsizes(time)/),float)
;ClO + BrO = BrCl
  data(0,:)=2*LOxLatMean(R_BrO_ClOa_col_x,j)
;ClO + BrO = Br + Cl
  data(1,:)=2*LOxLatMean(R_BrO_ClOb_col_x,j)
;JCl2O2 [Cl2O2]
  data(2,:)=2*LOxLatMean(J_Cl2O2_col_x,j)
;JHOBr [HOBr]
  data(3,:)=2*LOxLatMean(J_HOBr_col_x,j)
;  res@xyLineThicknessF = 2.0
;  delete(res@xyLineColors)
;  res@xyLineColors      = (/"blue", "red", "green", "black"/); change line color
  res@pmLegendDisplayMode = "NoCreate"
;  delete(res@xyExplicitLegendLabels)
;  res@xyExplicitLegendLabels=(/"2 x O~B~3~N~ + O", "2 x ClO + O", "2 x NO~B~2~N~ + O", "NO~D~3~N~ + hv", "HO2 + O", "HO2 + O3", "H + O3", "OH + O", "ClO + BrO -> BrCl", "ClO + BrO -> Br + Cl", "Cl2O2+hv", "BrO + O", "BrO + BrO", "HOBr+hv"/)
  plot(3) = gsn_csm_xy(wks, time, data, res)  

  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(plot)

;------Page 6--------    
  nPanel=4
  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)
  
  res@xyLineThicknessF = 4.0
  
  if (j.eq.-1) then
    resP@txString = "Global Mean Ox Production"
  else
    resP@txString = "Latitude: "+lat(j)+" Zonal Mean Ox Production"
  end if
  
  delete(res@xyLineColors)
  res@xyLineColors      = (/"blue","red"/)          ; change line color
  
  delete(data)
  data= new((/2,dimsizes(time)/),float)
    
  res@gsnRightString = "molec/cm~S~2~N~/s"
  delete(res@xyExplicitLegendLabels)
  res@xyExplicitLegendLabels=(/"Geoeng","Control"/)
    
  res@pmLegendDisplayMode = "Always"
  res@gsnLeftString = "O2 + hv -> 2O rate"
  data(0,:)=LOxLatMean(J_O2_b_col_x,j)
  data(1,:)=LOxLatMean(J_O2_b_col_c,j)
  plot(0) = gsn_csm_xy(wks, time, data, res)
    
  res@pmLegendDisplayMode = "NoCreate"
  res@gsnLeftString = "O2 + hv -> O + O(1D) rate"
  data(0,:)=LOxLatMean(J_O2_a_col_x,j)
  data(1,:)=LOxLatMean(J_O2_a_col_c,j)
  plot(1) = gsn_csm_xy(wks, time, data, res)

; see Brasseur and Solomon, p. 296
  res@gsnLeftString = "2 x CH3O2+NO"
  data(0,:)=2*LOxLatMean(R_CH3O2_NO_col_x,j)
  data(1,:)=2*LOxLatMean(R_CH3O2_NO_col_c,j)
  plot(2) = gsn_csm_xy(wks, time, data, res)
    
  res@gsnLeftString = "NO+HO2"
  data(0,:)=LOxLatMean(R_NO_HO2_col_x,j)
  data(1,:)=LOxLatMean(R_NO_HO2_col_c,j)
  plot(3) = gsn_csm_xy(wks, time, data, res)

  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(plot)

;------Page 7--------  
  
  nPanel=5
  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)
  
  if (j.eq.-1) then
    resP@txString = "Global Mean Ox Loss"
  else
    resP@txString = "Latitude: "+lat(j)+" Zonal Mean Ox Loss"
  end if

; NO + O3 --> NO2 + O2
; NO2 + O --> NO + O2
; ---------------------------
; NET: O3 + O --> 2 O2
  res@pmLegendDisplayMode = "Always"
  res@gsnLeftString = "2 x NO2+O rate"
  data(0,:)=2*LOxLatMean(R_NO2_O_col_x,j)
  data(1,:)=2*LOxLatMean(R_NO2_O_col_c,j)
  plot(0) = gsn_csm_xy(wks, time, data, res)

  res@pmLegendDisplayMode = "NoCreate"
  res@gsnLeftString = "2 x O3+O rate"
  data(0,:)=2*LOxLatMean(R_O3_O_col_x,j)
  data(1,:)=2*LOxLatMean(R_O3_O_col_c,j)
  plot(1) = gsn_csm_xy(wks, time, data, res)
  
; H + O2 + M --> HO2 + M
; HO2 + O --> OH + O2
; OH + O --> H + O2
; --------------------------
; NET: 2O --> O2  
  res@gsnLeftString = "OH+O rate"
  data(0,:)=LOxLatMean(R_OH_O_col_x,j)
  data(1,:)=LOxLatMean(R_OH_O_col_c,j)
  plot(2) = gsn_csm_xy(wks, time, data, res)
  
; OH + O3 --> HO2 + O2
; HO2 + O --> OH + O2
; ----------------------------
; NET: O3 + O --> 2 O2
  res@gsnLeftString = "HO2+O rate"
  data(0,:)=LOxLatMean(R_HO2_O_col_x,j)
  data(1,:)=LOxLatMean(R_HO2_O_col_c,j)
  plot(3) = gsn_csm_xy(wks, time, data, res)
    
  res@gsnLeftString = "H+O3 rate"
  data(0,:)=LOxLatMean(R_H_O3_col_x,j)
  data(1,:)=LOxLatMean(R_H_O3_col_c,j)
  plot(4) = gsn_csm_xy(wks, time, data, res)
  
  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(plot)

;------Page 8--------    
  nPanel=5
  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)
  
  res@pmLegendDisplayMode = "Always"
  res@gsnLeftString = "O2+O rate"
  data(0,:)=LOxLatMean(R_O2_O_col_x,j)
  data(1,:)=LOxLatMean(R_O2_O_col_c,j)
  plot(0) = gsn_csm_xy(wks, time, data, res)

; Cl + O3 --> ClO + O2
; ClO + O --> Cl + O2
; -------------------
; NET: O3 + O --> 2 O2
  res@pmLegendDisplayMode = "NoCreate"
  res@gsnLeftString = "2 x ClO+O rate"
  data(0,:)=2*LOxLatMean(R_ClO_O_col_x,j)
  data(1,:)=2*LOxLatMean(R_ClO_O_col_c,j)
  plot(1) = gsn_csm_xy(wks, time, data, res)


; ClONO2 + hv --> Cl + NO3
; NO3 + hv --> NO + O2
; NO + O3 --> NO2 + O2
; Cl + O3 --> ClO + O2
; ClO + NO2 + M --> ClONO2 + M
; ----------------------------
; NET: 2O3 --> 3O2
;
; BrO + NO2 + M --> BrONO2 + M
; BrONO2 + hv --> Br + NO3
; NO3 + hv --> NO + O2
; NO + O3 --> NO2 + O2
; Br + O3 --> BrO + O2
; --------------------
; NET: 2 O3 --> 3 O2        (
  res@gsnLeftString = "2 x NO3+hv rate"
  data(0,:)=2*LOxLatMean(J_NO3_b_col_x,j)
  data(1,:)=2*LOxLatMean(J_NO3_b_col_c,j)
  plot(2) = gsn_csm_xy(wks, time, data, res)

; OH + O -> H + O2
; H + O2 + M -> HO2 +M
; HO2 + O3 -> OH + 2O2
; --------------------
; NET: O + O3 -> 2 O2    
  res@gsnLeftString = "HO2+O3 rate"
  data(0,:)=LOxLatMean(R_HO2_O3_col_x,j)
  data(1,:)=LOxLatMean(R_HO2_O3_col_c,j)
  plot(3) = gsn_csm_xy(wks, time, data, res)
  
  res@gsnLeftString = "HOCl + hv rate"
  data(0,:)=LOxLatMean(J_HOCl_col_x,j)
  data(1,:)=LOxLatMean(J_HOCl_col_c,j)
  plot(4) = gsn_csm_xy(wks, time, data, res)
  
  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(plot)

;------Page 9--------    
  nPanel=6
  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)
  
  res@pmLegendDisplayMode = "Always"
  res@gsnLeftString = "2 x BrO + O rate"
  data(0,:)=2*LOxLatMean(R_BrO_O_col_x,j)
  data(1,:)=2*LOxLatMean(R_BrO_O_col_c,j)
  plot(0) = gsn_csm_xy(wks, time, data, res)
  
  res@pmLegendDisplayMode = "NoCreate"
  res@gsnLeftString = "2 x BrO + Bro rate"
  data(0,:)=2*LOxLatMean(R_BrO_BrO_col_x,j)
  data(1,:)=2*LOxLatMean(R_BrO_BrO_col_c,j)
  plot(1) = gsn_csm_xy(wks, time, data, res)
    
; ClO + BrO -> Cl + Br + O2
; Cl + O3 -> ClO + O2
; Br + O3 -> BrO + O2
; -------------------------
; NET: 2 O3 -> 3 O2
  res@gsnLeftString = "2 x ClO + BrO -> Br + Cl + O2 rate"
  data(0,:)=2*LOxLatMean(R_BrO_ClOa_col_x,j)
  data(1,:)=2*LOxLatMean(R_BrO_ClOa_col_c,j)
  plot(2) = gsn_csm_xy(wks, time, data, res)
    
; ClO + BrO -> BrCl + O2
; BrCl + hv -> Cl + Br
; Cl + O3 -> ClO + O2
; Br + O3 -> BrO + O2
; -------------------------
; NET: 2 O3 -> 3 O2
  res@gsnLeftString = "2 x ClO + BrO -> BrCl + O2 rate"
  data(0,:)=2*LOxLatMean(R_BrO_ClOb_col_x,j)
  data(1,:)=2*LOxLatMean(R_BrO_ClOb_col_c,j)
  plot(3) = gsn_csm_xy(wks, time, data, res)
    
; BrO + HO2 --> HOBr + O2
; HOBr + hv --> OH + Br
; OH + O3 --> HO2 + O2
; Br + O3 --> BrO + O2
; --------------------
; NET: 2 O3 --> 3 O2
  res@gsnLeftString = "2 x HOBr + hv rate"
  data(0,:)=2*LOxLatMean(J_HOBr_col_x,j)
  data(1,:)=2*LOxLatMean(J_HOBr_col_c,j)
  plot(4) = gsn_csm_xy(wks, time, data, res)

; ClO + ClO + M --> Cl2O2 + M
; Cl2O2 + hv --> Cl + ClO2
; ClOO + M --> Cl + O2 + M
; 2 (Cl + O3 --> ClO + O2)
; ----------------------------
; NET: 2 O3 --> 3 O2    
  res@gsnLeftString = "2 x Cl2O2 + hv rate"
  data(0,:)=2*LOxLatMean(J_Cl2O2_col_x,j)
  data(1,:)=2*LOxLatMean(J_Cl2O2_col_c,j)
  plot(5) = gsn_csm_xy(wks, time, data, res)
  
  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(wks) 
    
end

; Main
; example: plotLox("experiment","control",j)
; j = latitude index, or -1 for global mean

;plotLOx("geoengn4","carmasadn3")
plotLOx("geo2","nongeo2",-1)
;plotLOx("geo2","nongeo2",35) ; 50N
;plotLOx("geo2","nongeo2",44) ; 86N
;plotLOx("geo2","nongeo2",1) ; 86S
;plotLOx("geo2","nongeo2",10) ; 50N
