load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

season="DJF"
;season="JJA"

infile=addfile("Ujet_3panel_std_wcm.nc","r")
Udiff_wcm=infile->Udiff
alpha_wcm=infile->alpha
Uavg_wcm60=infile->Uavg
delete(infile)

infile=addfile("Ujet_3panel_std_c1d.nc","r")
Udiff_c1d=infile->Udiff
alpha_c1d=infile->alpha
Uavg_c1d60=infile->Uavg
delete(infile)

infile=addfile("Ujet_3panel_std_wset.nc","r")
Udiff_wset=infile->Udiff
alpha_wset=infile->alpha
Uavg_wset60=infile->Uavg
delete(infile)

  res = True
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@cnLinesOn             = False
  res@cnFillOn              = True                  ; turn on color fill
  res@gsnSpreadColors     = True ; spread out color table
  res@cnLevelSelectionMode  = "ManualLevels" 
; Contours by 10 up to 200:
  res@cnLevelSpacingF       = 1.0                  ; contour interval
  res@cnMinLevelValF        = -10.0                  ; contour interval
  res@cnMaxLevelValF        = 10.0                 ; contour interval
;  res@lbLabelAngleF         = 45
  res@lbLabelStride         = 2
  
  res@tmYLMode      = "Explicit" ; Define own tick mark labels.
;  res@tmYLLabelFont = 25         ; Change the font.
  res@tmYLValues    = (/1000,300,100,30,10,3,1,0.3,0.1, 0.01, 0.001, 0.0001/)
  res@tmYLLabels    = (/"1000","300","100","30","10","3","1","0.3","0.1", "0.01", "0.001", "0.0001"/)
  res@tmYLMinorOn   = True        ; No minor tick marks
  res@lbLabelBarOn  = False
  res@cnLineLabelsOn       = False       ; do not draw contour labels
  
  res2 = True
  res2@gsnDraw              = False
  res2@gsnFrame             = False
  res2@cnLevelSelectionMode = "ExplicitLevels" ; set manual contour levels
  res2@cnLevels             = (/ -5.0, 95.0, 105.0 /)
  res2@cnInfoLabelOn        = False       ; turn off info label
  res2@cnLinesOn            = False       ; do not draw contour lines
  res2@cnFillOn             = False
  res2@cnLineLabelsOn       = False       ; do not draw contour labels
;  res2@cnFillScaleF         = 0.6         ; add extra density
;  res2@cnMonoFillColor      = False
;  res2@cnFillColors         = (/ Foreground, Black /)
;  res2@cnMonoFillPattern    = False
;  res2@cnFillPatterns       = (/ 4, 17 /)   ; shade areas less than 95% confidence
  
  res3 = True
  res3@gsnDraw               = False
  res3@gsnFrame              = False
  res3@cnLinesOn             = True
  res3@cnFillOn              = False
  res3@cnLevelSelectionMode  = "ExplicitLevels" 
  cnLevels = ispan(-30, 45, 5)
  cnLineDashPatterns         = where(cnLevels.lt.0, 1, 0)
  cnLineThicknesses         = where(cnLevels.eq.0, 2, 1)
; Contours -30 to 45 by 5:
  res3@cnLevels              = cnLevels
  res3@cnMonoLineDashPattern = False
  res3@cnLineDashPatterns    = cnLineDashPatterns
  res3@cnMonoLineThickness   = False
  res3@cnLineThicknesses     = cnLineThicknesses
  
  imagefile="$IMAGE_PATH/3case_Ujet_"+season
  wks = gsn_open_wks("ps",imagefile)
  gsn_define_colormap(wks,"blwhre")                ; choose colormap
  
  minP = 10
  maxP = 1013.
  
; Color filled contours of wind shifts
  plot1 = new(3,graphic)  
  plot1(0) = gsn_csm_pres_hgt(wks, Udiff_wcm({minP:maxP}, {-90.0:90.0}), res)
  plot1(1) = gsn_csm_pres_hgt(wks, Udiff_c1d({minP:maxP}, {-90.0:90.0}), res)
  plot1(2) = gsn_csm_pres_hgt(wks, Udiff_wset({minP:maxP}, {-90.0:90.0}), res)

; Add significance stippling
  plot2 = new(3,graphic)  
  plot2(0) = gsn_csm_pres_hgt(wks, alpha_wcm({minP:maxP}, {-90.0:90.0}), res2)
  plot2(1) = gsn_csm_pres_hgt(wks, alpha_c1d({minP:maxP}, {-90.0:90.0}), res2)
  plot2(2) = gsn_csm_pres_hgt(wks, alpha_wset({minP:maxP}, {-90.0:90.0}), res2)
  opt = True
  opt@gsnShadeFillType = "pattern"      ; color is the default
  opt@gsnShadeLow = 17                 ; Fill pattern 17 is stippled  
  ; see http://www.ncl.ucar.edu/Document/Graphics/Images/fillpatterns.png
  plot2(0) = gsn_contour_shade (plot2(0), 95., 105., opt)
  plot2(1) = gsn_contour_shade (plot2(1), 95., 105., opt)
  plot2(2) = gsn_contour_shade (plot2(2), 95., 105., opt)

; Add line contours of winds for 1960-70
  plot3 = new(3,graphic)  
  plot3(0) = gsn_csm_pres_hgt(wks, Uavg_wcm60({minP:maxP}, {-90.0:90.0}), res3)
  plot3(1) = gsn_csm_pres_hgt(wks, Uavg_c1d60({minP:maxP}, {-90.0:90.0}), res3)
  plot3(2) = gsn_csm_pres_hgt(wks, Uavg_wset60({minP:maxP}, {-90.0:90.0}), res3)

  overlay(plot1(0), plot2(0))
  overlay(plot1(0), plot3(0))
  overlay(plot1(1), plot2(1))
  overlay(plot1(1), plot3(1))
  overlay(plot1(2), plot2(2))
  overlay(plot1(2), plot3(2))

  resP                            = True
  resP@gsnMaximize                = True
  resP@gsnPanelXWhiteSpacePercent = 3
  resP@gsnPaperOrientation        = "landscape"
  resP@gsnPanelLabelBar           = True
  
  gsn_panel(wks,plot1,(/1,3/),resP)
;  draw(plot1)
  frame(wks)
  delete(wks)
  
  muttcmd="mutt mikejmills@mac.com -a "+imagefile+".ps -s "+imagefile+" < /dev/null"
  print(muttcmd)
  print("open "+imagefile+".ps")
;  mutter=systemfunc(muttcmd)
;  print(mutter)
