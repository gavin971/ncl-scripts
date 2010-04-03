load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("plotvgO3Diff")
procedure plotvgO3Diff(O3_col_x[*][*]:float, O3_col_c[*][*]:float, \
   O3_colSD_x[*][*]:float, O3_colSD_c[*][*]:float)
   
local O3_col_diff,wks,res,lev
begin 

  lat=O3_col_x&lat
  nLat=(/dimsizes(lat)/)
  
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
  
  iflag= False               ; population variance similar
  prob = ttest(O3_col_x, O3_colSD_x, 30, O3_col_c, O3_colSD_c, 30, iflag, False)
  prob!0 = "time"
  prob!1 = "lat"
  prob&lat = lat
  
  sig_level = 0.025  ; significance level (95%)
  
  fmt = "%5.2f"
  monthnames = (/ "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"/)
;  monthnames=(/"J","F","M","A","M","J","J","A","S","O","N","D"/)
  
  print((/lat(0)/))
  do i = 1,ntimes -2
    print(monthnames(i-1)+"  "+sprintf(fmt, O3_col_x(i,0))+"  "+sprintf(fmt, O3_col_c(i,0))+ \
       "  "+sprintf(fmt, O3_col_diff(i,0))+"  "+sprintf(fmt, O3_col_dpct(i,0))+"  "+\
       sprintf(fmt, prob(i,0)))
  end do
  
  print((/lat(nLat-1)/))
  do i = 1,ntimes -2
    print(monthnames(i-1)+"  "+sprintf(fmt, O3_col_x(i,nLat-1))+"  "+sprintf(fmt, O3_col_c(i,nLat-1))+ \
       "  "+sprintf(fmt, O3_col_diff(i,nLat-1))+"  "+sprintf(fmt, O3_col_dpct(i,nLat-1))+"  "+\
       sprintf(fmt, prob(i,nLat-1)))
  end do
  
;  print("probability:")
;  do i = 1,ntimes -2
;    print(monthnames(i-1)+"  "+sprintf(fmt, (/prob(i,:)/)))
;  end do
  
;  print("t-values:")
;  fmt = "%5.2f"
;  do i = 1,ntimes -2
;    print(monthnames(i-1)+"  "+sprintf(fmt, prob(1,i,:)))
;  end do
  
  if (.not.(isfilepresent("images/vg"))) then
    a = systemfunc("mkdir -p images/vg")
  end if  
  wks = gsn_open_wks("pdf", "images/vg/vgO3ColumnChange")
    
  resP                      = True
  resP@gsnMaximize          = True
  resP@gsnPanelXWhiteSpacePercent = 5
  resP@gsnPaperOrientation = "landscape"
  
  res                     = True          ; plot mods desired
  res@gsnPaperOrientation = "landscape"

  res@cnFillOn            = True          ; turn on color
;  res@lbLabelBarOn        = False         ; turn off color labelbar
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.7
  res@gsnLeftStringFontHeightF = 24
  res@gsnRightStringFontHeightF = 24
  res@tiYAxisString       = "Latitude"
  res@cnLineLabelPlacementMode= "Computed"
  res@cnLineLabelDensityF = 5
  res@tmXBLabelFontHeightF =0.025
  res@tmYLLabelFontHeightF =0.03
  
  res@tmXBLabels = monthnames
  res@tmXBMode           = "Explicit"         ; Define own tick mark labels.
  res@tmXBValues         = O3_col_c&time(1:12)
  
  res@tmXBMinorValues    = O3_col_c&time(0:12)+0.5
  res@tmXBMinorOutwardLengthF = 0.02    ; Increase length and size of
  res@tmXBMinorLengthF        = 0.02    ; minor tickmarks
  res@tmXBMinorThicknessF     = 2.0
;  res@tmXBMajorThicknessF     = 2.0
;  res@tmXBMajorOutwardLengthF = 0.01
  res@tmXBMajorLengthF        = 0.00
;===================================================================  
  res@trXMinF = 0.5   
  res@trXMaxF = 12.5  
; Define our own contour levels.  
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnMinLevelValF=150
  res@cnMaxLevelValF=450
  res@cnLevelSpacingF     = 30.
  res@gsnRightString      = "O~B~3~N~ Column(DU)"
  res@lbLabelAngleF       = 60
  res@gsnLeftString       = "Control run"
  
  gsn_define_colormap(wks,"rainbow")
  plot = gsn_csm_lat_time(wks,O3_col_c(lat|:, time|:),res)      ; contour the variable
  
  res@gsnLeftString       = "Soot run"
  plot = gsn_csm_lat_time(wks,O3_col_x(lat|:, time|:),res)      ; contour the variable
  
;===================================================================  

  nPlots = 1
  plot = new((/ nPlots /), graphic)
  plot2 = new((/ nPlots /), graphic)
    
  res@gsnDraw               = False       ; Do not draw plot
  res@gsnFrame              = False       ; Do not advance frame
  
;===================================================================  
  
  res2 = True                            ; res2 probability plots
  res2@gsnDraw             = False       ; Do not draw plot
  res2@gsnFrame            = False       ; Do not advance frame
  res2@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res2@cnMinLevelValF      = 0.00        ; set min contour level
  res2@cnMaxLevelValF      = 1.00        ; set max contour level
  res2@cnLevelSpacingF     = sig_level   ; set contour spacing
  res2@cnInfoLabelOn       = False       ; turn off info label
  res2@cnLinesOn           = False       ; do not draw contour lines
  res2@cnLineLabelsOn      = False       ; do not draw contour labels
  res2@cnFillScaleF        = 0.6         ; add extra density

;===================================================================  

  gsn_define_colormap(wks,"BlRe")
  res@cnLevelSelectionMode = "ExplicitLevels" 
  res@cnLevels          = (/ -4, -3, -2, -1,-0.5, 0, 1, 2, 5, 10, 20/) 
  res@cnLineThicknesses = (/ 1 , 1 , 1 , 1 ,   1, 4, 1, 1, 1, 1,   1 /)
  res@cnMonoLineThickness = False
  res@cnLevelSpacingF     = 30.
  res@gsnLeftString       = "Soot run - Control run"
  
  plot = gsn_csm_lat_time(wks, O3_col_diff(lat|:, time|:),res)      ; contour the variable
  opt = True
  opt@gsnShadeFillType = "pattern"      ; pattern fill
  opt@gsnShadeHigh = 17                  ; pattern # to use 
  plot2 = gsn_csm_lat_time(wks, prob(lat|:, time|:), res2) 
  plot2   = gsn_contour_shade(plot2, -999, sig_level, opt)  ;
  overlay (plot, plot2)
  draw(plot)
  frame(wks)
  
  res@cnLevelSelectionMode = "ExplicitLevels" 
  delete(res@cnLevels)
  res@cnLevels          = (/ -4, -3, -2 , -1,-0.5,-.25, 0, 0.5, 1, 2, 3, 5, 10/) 
  delete(res@cnLineThicknesses)
  res@cnLineThicknesses = (/  1 , 1 , 1 , 1 ,  1 ,   1, 4,   1, 1, 1, 1, 1, 1 /)
  res@cnMonoLineThickness = False
  res@gsnLeftString       = "(Soot - Control)/Control"
  res@gsnRightString      = "O~B~3~N~ Column Change(%)"
  plot = gsn_csm_lat_time(wks,O3_col_dpct(lat|:, time|:),res)      ; contour the variable
  overlay (plot, plot2)
  draw(plot)
  frame(wks)
end


undef("plotvgO3Change")
procedure plotvgO3Change(rootPath:string,historyRun:string, controlRun:string)
local T_x,T_c

begin
  print((/historyRun/))
;  saveMeanVars(historyRun,True)
;  saveMeanLOxFiles(historyRun)
  historyPattern = rootPath+"/"+historyRun+"/h0zm/"+historyRun+".cam2.h0zm.[01][0-9].nc"
  historyNames = systemfunc("ls " + historyPattern)
;  print((/historyNames/))
  historySDPattern = rootPath+"/"+historyRun+"/h0zm/"+historyRun+".cam2.h0zmsd.[01][0-9].nc"
  historySDNames = systemfunc("ls " + historySDPattern)
  
  print((/controlRun/))
  controlPattern = rootPath+"/"+controlRun+"/h0zm/"+controlRun+".cam2.h0zm.[01][0-9].nc"
  controlNames = systemfunc("ls " + controlPattern)
;  print((/controlNames/))
  controlSDPattern = rootPath+"/"+controlRun+"/h0zm/"+controlRun+".cam2.h0zmsd.[01][0-9].nc"
  controlSDNames = systemfunc("ls " + controlSDPattern)
  
  a=addfile(historyNames(0),"r")
  O3_col_x=new((/14,dimsizes(a->lat)/),float)
  O3_col_x@units="Dobson Units"
  O3_col_x!0 = "time"
  O3_col_x!1 = "lat"
  O3_col_x&lat = a->lat
  print(dimsizes(O3_col_x))
  O3_col_c=O3_col_x
  O3_colSD_x=O3_col_x
  O3_colSD_c=O3_col_x
  
  time = new(14,integer)
  
  do i = 1,12
    time(i)=i
    a=addfile(historyNames(i-1),"r")
    b=addfile(controlNames(i-1),"r")
    c=addfile(historySDNames(i-1),"r")
    d=addfile(controlSDNames(i-1),"r")
    O3_col_x(i,:)=a->O3_col
    O3_col_c(i,:)=b->O3_col
    O3_colSD_x(i,:)=c->O3_col
    O3_colSD_c(i,:)=d->O3_col
  end do
  time(13)=13
  time(0)=0
  
  O3_col_x(0,:)=O3_col_x(12,:)
  O3_col_c(0,:)=O3_col_c(12,:)
  O3_colSD_x(0,:)=O3_col_x(12,:)
  O3_colSD_c(0,:)=O3_col_c(12,:)
  
  O3_col_x(13,:)=O3_col_x(1,:)
  O3_col_c(13,:)=O3_col_c(1,:)
  O3_colSD_x(13,:)=O3_col_x(1,:)
  O3_colSD_c(13,:)=O3_col_c(1,:)
  
  O3_col_x&time=time
  O3_col_c&time=time
  O3_colSD_x&time=time
  O3_colSD_c&time=time

  print("plotting ozone column...")
  plotvgO3Diff(O3_col_x, O3_col_c, O3_colSD_x, O3_colSD_c)
  
end

historyRun="vgsoot08"
controlRun="vgcontrol02"
;rootPath="/data5/mmills"
;rootPath="/Volumes/Data/Models/waccm/carma-soot/run"
rootPath="/Data/Models/waccm/carma-soot/run"

plotvgO3Change(rootPath,historyRun,controlRun)
