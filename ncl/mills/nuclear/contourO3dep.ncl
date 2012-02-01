load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("plotO3Diff")
procedure plotO3Diff(historyRun:string, O3_col_x[*][*]:float, O3_col_c[*][*]:float) 
local O3_col_diff,wks,res,lev
begin 

  lat=O3_col_x&lat
  
  xdim=dimsizes(O3_col_x)
  cdim=dimsizes(O3_col_c)
  print("xdim="+xdim+" cdim="+cdim)
  ntimes=min((/xdim(0),cdim(0)/))
  
  time=O3_col_x&time
  time=time/365.  
  
  O3_col_diff=O3_col_x(:ntimes-1,:)-O3_col_c(:ntimes-1,:)
  O3_col_diff!0 = "time"
  O3_col_diff!1 = "lat"
  O3_col_diff&lat = lat
  O3_col_diff&time = time

  O3_col_dpct=100.*O3_col_diff(:ntimes-1,:)/O3_col_c(:ntimes-1,:)
  O3_col_dpct!0 = "time"
  O3_col_dpct!1 = "lat"
  O3_col_dpct&lat = lat
  O3_col_dpct&time = time
  
  
  pdffile = "$IMAGE_PATH/ContourO3dep."+historyRun
  wks = gsn_open_wks("pdf", pdffile)
  gsn_define_colormap(wks,"posneg_1")
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
  res@tiXAxisString       = ""
  res@tiYAxisString       = "Latitude"
  res@cnLineLabelPlacementMode= "Computed"
  res@cnLineLabelDensityF = 5
  res@tmXBLabelFontHeightF =0.015
  res@tmYLLabelFontHeightF =0.03
 
  res@tmXBMaxTicks      = 10
  res@tmXBMinorPerMajor = 3
  res@tmXBAutoPrecision = False
  res@tmXBPrecision = 2

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
  res@cnLevels          = (/-70, -60,-50,-40,-30,-20,-10,-5,-2,0,2,5/)
  delete(res@cnLineThicknesses)
  res@cnLineThicknesses = (/ 1,  1 , 1 , 1 , 1 , 1 , 1 , 1, 1,4, 1, 1/)
  res@cnMonoLineThickness = False
  res@gsnLeftString       = ""
  res@gsnRightString      = ""
  res@gsnCenterString      = "O3 Column Change (%)"
  plot = gsn_csm_lat_time(wks,O3_col_dpct(lat|:, time|:),res)      ; contour the variable
  delete(plot)
  
  ostype = systemfunc("echo $OSTYPE")
  if (ostype.eq."darwin") then
    opencmd="open "+pdffile+".pdf"
    print(opencmd)
    open=systemfunc(opencmd)
    print(open)
  end if
end

undef("plotO3Change")
procedure plotO3Change(rootPath:string, historyRun:string, controlRun:string)
local T_x,T_c

begin
  print((/historyRun/))
  historyPath2=rootPath+historyRun+"/"
  analysisPath=historyPath2+"proc/"
  historyPattern = analysisPath+historyRun+".O3Col.h0zm.nc"
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  print((/historyNames/))
  O3_col_x=addfiles_GetVar(history, historyNames, "O3_col")
  print(dimsizes(O3_col_x))
 
;  T_x=File1->T
;  Soot=File1->SOOT01

  print((/controlRun/))
  controlPath2=rootPath+controlRun+"/"
  analysisPath=controlPath2+"proc/"
  controlPattern = analysisPath+controlRun+".O3Col.h0zm.nc"
  controlNames = systemfunc("ls " + controlPattern)
  control = addfiles(controlNames, "r")
  ListSetType(control, "cat")
  print((/controlNames/))
  O3_col_c=addfiles_GetVar(control, controlNames, "O3_col")
  print(dimsizes(O3_col_c))
  
  print("plotting ozone column...")
  plotO3Diff(historyRun, O3_col_x, O3_col_c)
    
end

host   = systemfunc("hostname")
if (host.eq."modeling1") then
  rootPath="/data4/mmills/"
else
  rootPath="$MODEL_PATH/ccsm/case/"
end if

historyRun="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"
controlRun="b40.rcp4_5.2deg.wcm.carma.bc0tg.avg"
plotO3Change(rootPath, historyRun, controlRun)
