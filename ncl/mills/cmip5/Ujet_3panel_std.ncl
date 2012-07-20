load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

host=systemfunc("hostname")
;print(host)

season="DJF"
;season="JJA"

; mid 21stC shift vs late 20thC
;rootdir="/data/waccm"
;case1="b40.1955-2005.2deg.wcm.avg"
;case2="b40.rcp4_5.2deg.wcm.avg"
;dir1="h0zm.seasavg.1986-2005"
;dir2="h0zm.seasavg.2032-2051"

;=========
; Late 20thC WACCM shift
;rootdir="/data/waccm"
;case1="b40.1955-2005.2deg.wcm.avg"
;case2=case1
;dir1="h0zm.seasavg.1960-1979"
;dir2="h0zm.seasavg.1986-2005"
rootdir = "/Volumes/MillsPassport/cesm/case/"
cases_wcm="b40.1955-2005.2deg.wcm.00[234]"
infiles_wcm60 = systemfunc("ls "+rootdir+cases_wcm+"/proc/*.cam2.h0zm.U."+season+".1960-1979.nc")
history_wcm60 = addfiles(infiles_wcm60, "r")
ListSetType(history_wcm60, "cat")
U_wcm60 = addfiles_GetVar(history_wcm60, infiles_wcm60, "U")
Uavg_wcm60=dim_avg_Wrap(U_wcm60(lev|:,lat|:,time|:))
Uvar_wcm60=dim_variance_Wrap(U_wcm60(lev|:,lat|:,time|:))

lev_wcm60=Uavg_wcm60&lev
lev_wcm60@units="hPa"
lev_wcm60@long_name="Pressure"

Uavg_wcm60&lev=lev_wcm60
Uvar_wcm60&lev=lev_wcm60

infiles_wcm86 = systemfunc("ls "+rootdir+cases_wcm+"/proc/*.cam2.h0zm.U."+season+".1986-2005.nc")
history_wcm86 = addfiles(infiles_wcm86, "r")
ListSetType(history_wcm86, "cat")
U_wcm86 = addfiles_GetVar(history_wcm86, infiles_wcm86, "U")
Uavg_wcm86=dim_avg_Wrap(U_wcm86(lev|:,lat|:,time|:))
Uvar_wcm86=dim_variance_Wrap(U_wcm86(lev|:,lat|:,time|:))
Uavg_wcm86&lev=lev_wcm60
Uvar_wcm86&lev=lev_wcm60

s_wcm60   = dimsizes (U_wcm60(time|:,lev|0,lat|0)) ; time sample sizes
s_wcm86   = dimsizes (U_wcm86(time|:,lev|0,lat|0))

iflag= False               ; population variance similar
siglvl = 0.05
prob_wcm = ttest(Uavg_wcm60,Uvar_wcm60,s_wcm60, Uavg_wcm86,Uvar_wcm86,s_wcm86, iflag, False) 
alpha_wcm = Uavg_wcm60 ; get metadata
alpha_wcm = (/100.*(1. - prob_wcm)/) ; A significance of 0.05 returned by ttest would yield 95% for alpha

Udiff_wcm=Uavg_wcm86 ; get metadata
Udiff_wcm=Udiff_wcm-Uavg_wcm60


;=========
; Late 20thC CAM4 1-deg shift
;rootdir="/data/ccsm4"
rootdir = "/Volumes/MillsPassport/cesm/case/"
cases_c1d="b40.20th.track1.1deg.00[5-9]"
infiles_c1d60 = systemfunc("ls "+rootdir+cases_c1d+"/proc/*.cam2.h0zm.U."+season+".1960-1979.nc")
history_c1d60 = addfiles(infiles_c1d60, "r")
ListSetType(history_c1d60, "cat")
U_c1d60 = addfiles_GetVar(history_c1d60, infiles_c1d60, "U")
Uavg_c1d60=dim_avg_Wrap(U_c1d60(lev|:,lat|:,time|:))
Uvar_c1d60=dim_variance_Wrap(U_c1d60(lev|:,lat|:,time|:))

lev_c1d60=Uavg_c1d60&lev
lev_c1d60@units="hPa"
lev_c1d60@long_name="Pressure"

Uavg_c1d60&lev=lev_c1d60
Uvar_c1d60&lev=lev_c1d60

infiles_c1d86 = systemfunc("ls "+rootdir+cases_c1d+"/proc/*.cam2.h0zm.U."+season+".1986-2005.nc")
history_c1d86 = addfiles(infiles_c1d86, "r")
ListSetType(history_c1d86, "cat")
U_c1d86 = addfiles_GetVar(history_c1d86, infiles_c1d86, "U")
Uavg_c1d86=dim_avg_Wrap(U_c1d86(lev|:,lat|:,time|:))
Uvar_c1d86=dim_variance_Wrap(U_c1d86(lev|:,lat|:,time|:))
Uavg_c1d86&lev=lev_c1d60
Uvar_c1d86&lev=lev_c1d60

s_c1d60   = dimsizes (U_c1d60(time|:,lev|0,lat|0)) ; time sample sizes
s_c1d86   = dimsizes (U_c1d86(time|:,lev|0,lat|0))

iflag= False               ; population variance similar
siglvl = 0.05
prob_c1d = ttest(Uavg_c1d60,Uvar_c1d60,s_c1d60, Uavg_c1d86,Uvar_c1d86,s_c1d86, iflag, False) 
alpha_c1d = Uavg_c1d60 ; get metadata
alpha_c1d = (/100.*(1. - prob_c1d)/) ; A significance of 0.05 returned by ttest would yield 95% for alpha

Udiff_c1d=Uavg_c1d86 ; get metadata
Udiff_c1d=Udiff_c1d-Uavg_c1d60

;=========
; Late 20thC wset shift
rootdir = "/Volumes/MillsPassport/cesm/case/"
cases_wset="b40.1955-2005.2deg.wset.00[456]"
infiles_wset60 = systemfunc("ls "+rootdir+cases_wset+"/proc/*.cam2.h0zm.U."+season+".1960-1979.nc")
history_wset60 = addfiles(infiles_wset60, "r")
ListSetType(history_wset60, "cat")
U_wset60 = addfiles_GetVar(history_wset60, infiles_wset60, "U")
Uavg_wset60=dim_avg_Wrap(U_wset60(lev|:,lat|:,time|:))
Uvar_wset60=dim_variance_Wrap(U_wset60(lev|:,lat|:,time|:))

lev_wset60=Uavg_wset60&lev
lev_wset60@units="hPa"
lev_wset60@long_name="Pressure"

Uavg_wset60&lev=lev_wset60
Uvar_wset60&lev=lev_wset60

infiles_wset86 = systemfunc("ls "+rootdir+cases_wset+"/proc/*.cam2.h0zm.U."+season+".1986-2005.nc")
history_wset86 = addfiles(infiles_wset86, "r")
ListSetType(history_wset86, "cat")
U_wset86 = addfiles_GetVar(history_wset86, infiles_wset86, "U")
Uavg_wset86=dim_avg_Wrap(U_wset86(lev|:,lat|:,time|:))
Uvar_wset86=dim_variance_Wrap(U_wset86(lev|:,lat|:,time|:))
Uavg_wset86&lev=lev_wset60
Uvar_wset86&lev=lev_wset60

s_wset60   = dimsizes (U_wset60(time|:,lev|0,lat|0)) ; time sample sizes
s_wset86   = dimsizes (U_wset86(time|:,lev|0,lat|0))

iflag= False               ; population variance similar
siglvl = 0.05
prob_wset = ttest(Uavg_wset60,Uvar_wset60,s_wset60, Uavg_wset86,Uvar_wset86,s_wset86, iflag, False) 
alpha_wset = Uavg_wset60 ; get metadata
alpha_wset = (/100.*(1. - prob_wset)/) ; A significance of 0.05 returned by ttest would yield 95% for alpha

Udiff_wset=Uavg_wset86 ; get metadata
Udiff_wset=Udiff_wset-Uavg_wset60



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
