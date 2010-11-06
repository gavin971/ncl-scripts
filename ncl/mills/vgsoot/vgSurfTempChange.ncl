load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("plotvgSurfTempDiff")
procedure plotvgSurfTempDiff(TS_x[*][*]:float, TS_c[*][*]:float, \
   TSSD_x[*][*]:float, TSSD_c[*][*]:float)
   
local TS_diff,wks,res,lev
begin 

  lat=TS_x&lat
  nLat=(/dimsizes(lat)/)
  
  xdim=dimsizes(TS_x)
  cdim=dimsizes(TS_c)
  print("xdim="+xdim+" cdim="+cdim)
  ntimes=min((/xdim(0),cdim(0)/))
  
  TS_diff=TS_x(:ntimes-1,:)-TS_c(:ntimes-1,:)
  TS_diff!0 = "time"
  TS_diff!1 = "lat"
  TS_diff&lat = lat
  
  iflag= False               ; population variance similar
  prob = ttest(TS_x, TSSD_x^2.0, 30, TS_c, TSSD_c^2.0, 30, iflag, False)
  prob!0 = "time"
  prob!1 = "lat"
  prob&lat = lat

  prob(0,:)=prob(12,:)
  prob(13,:)=prob(1,:)
  
;  sig_level = 0.025  ; significance level (95%)
  sig_level = 0.05  ; significance level (95%)
  
  fmt = "%5.2f"
  monthnames = (/ "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"/)
;  monthnames=(/"J","F","M","A","M","J","J","A","S","O","N","D"/)
  
  print((/lat(0)/))
  do i = 1,ntimes -2
    print(monthnames(i-1)+"  "+sprintf(fmt, TS_x(i,0))+"  "+sprintf(fmt, TS_c(i,0))+ \
       "  "+sprintf(fmt, TS_diff(i,0))+"  "+sprintf(fmt, prob(i,0)))
  end do
  
  print((/lat(nLat-1)/))
  do i = 1,ntimes -2
    print(monthnames(i-1)+"  "+sprintf(fmt, TS_x(i,nLat-1))+"  "+sprintf(fmt, TS_c(i,nLat-1))+ \
       "  "+sprintf(fmt, TS_diff(i,nLat-1))+"  "+sprintf(fmt, prob(i,nLat-1)))
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
  
  imgPath=getenv ("IMG_PATH")
  if (.not.(isfilepresent(imgPath+"/bc"))) then
    a = systemfunc("mkdir -p "+imgPath+"/bc")
  end if
  filename=imgPath+"/bc/bcSurfTempChange"
  wks = gsn_open_wks("ps", filename)
    
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
  res@tmXBValues         = TS_c&time(1:12)
  
  res@tmXBMinorValues    = TS_c&time(0:12)+0.5
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
  res@cnMinLevelValF=220
  res@cnMaxLevelValF=300
  res@cnLevelSpacingF     = 5.
  res@gsnRightString      = "Surface Temperature (K)"
  res@lbLabelAngleF       = 60
  res@gsnLeftString       = "Control run"
  
  gsn_define_colormap(wks,"rainbow")
  plot = gsn_csm_lat_time(wks,TS_c(lat|:, time|:),res)      ; contour the variable
  
  res@gsnLeftString       = "Soot run"
  plot = gsn_csm_lat_time(wks,TS_x(lat|:, time|:),res)      ; contour the variable
  
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
  res@cnLevels          = (/ -2, -1.5, -1, -0.5,-0.25, 0,  0.1, 0.25, 0.5, 0.75, 1.0 /) 
  res@cnLineThicknesses = (/ 1 ,   1 , 1 ,   1 ,   1,  4,    1,   1,    1, 1,   1 /)
  res@cnMonoLineThickness = False

  res@gsnLeftString       = "Soot run - Control run"
  
  plot = gsn_csm_lat_time(wks, TS_diff(lat|:, time|:),res)      ; contour the variable
  opt = True
  opt@gsnShadeFillType = "pattern"      ; pattern fill
  opt@gsnShadeHigh = 17                  ; pattern # to use 
  plot2 = gsn_csm_lat_time(wks, prob(lat|:, time|:), res2) 
  plot2   = gsn_contour_shade(plot2, -999, sig_level, opt)  ;
  overlay (plot, plot2)
  draw(plot)
  frame(wks)
  
  print("ps2pdf " + filename + ".ps" + " " + filename + ".pdf")
  pstest=systemfunc("ps2pdf " + filename + ".ps" + " " + filename + ".pdf")
  print(pstest)

  print("mutt -a " + filename + ".pdf mikejmills@mac.com < /dev/null")
  mutter=systemfunc("mutt -a " + filename + ".pdf mikejmills@mac.com < /dev/null")
  print(mutter)
  
end


undef("plotvgSurfTempChange")
procedure plotvgSurfTempChange(rootPath:string,historyRun:string, controlRun:string)
local TS_x,TS_c

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
  TS_x=new((/14,dimsizes(a->lat)/),float)
  TS_x@units="K"
  TS_x!0 = "time"
  TS_x!1 = "lat"
  TS_x&lat = a->lat
  print(dimsizes(TS_x))
  TS_c=TS_x
  TSSD_x=TS_x
  TSSD_c=TS_x
  
  time = new(14,integer)
  
  do i = 1,12
    time(i)=i
    a=addfile(historyNames(i-1),"r")
    b=addfile(controlNames(i-1),"r")
    c=addfile(historySDNames(i-1),"r")
    d=addfile(controlSDNames(i-1),"r")
    TS_x(i,:)=a->TS
    TS_c(i,:)=b->TS
    TSSD_x(i,:)=c->TS
    TSSD_c(i,:)=d->TS
  end do
  time(13)=13
  time(0)=0
  
  TS_x(0,:)=TS_x(12,:)
  TS_c(0,:)=TS_c(12,:)
  TSSD_x(0,:)=TS_x(12,:)
  TSSD_c(0,:)=TS_c(12,:)
  
  TS_x(13,:)=TS_x(1,:)
  TS_c(13,:)=TS_c(1,:)
  TSSD_x(13,:)=TS_x(1,:)
  TSSD_c(13,:)=TS_c(1,:)
  
  TS_x&time=time
  TS_c&time=time
  TSSD_x&time=time
  TSSD_c&time=time

  print("plotting ozone column...")
  plotvgSurfTempDiff(TS_x, TS_c, TSSD_x, TSSD_c)
  
end

historyRun="bcrckt03"
controlRun="bccntl02"
;rootPath="/data5/mmills"
;rootPath="/Volumes/Data/Models/waccm/carma-soot/run"
;rootPath="/Data/Models/waccm/run/carma-soot"
rootPath="/data5/mmills/wa319/carma-soot"

plotvgSurfTempChange(rootPath,historyRun,controlRun)
