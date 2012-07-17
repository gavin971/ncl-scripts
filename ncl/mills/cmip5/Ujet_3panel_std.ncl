load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

; mid 21stC shift vs late 20thC
;rootdir="/data/waccm"
;case1="b40.1955-2005.2deg.wcm.avg"
;case2="b40.rcp4_5.2deg.wcm.avg"
;dir1="h0zm.seasavg.1986-2005"
;dir2="h0zm.seasavg.2032-2051"

; Late 20thC WACCM shift
;rootdir="/data/waccm"
;case1="b40.1955-2005.2deg.wcm.avg"
;case2=case1
;dir1="h0zm.seasavg.1960-1979"
;dir2="h0zm.seasavg.1986-2005"

; Late 20thC WACCM shift
rootdir="/data/ccsm4"
case1="b40.1955-2005.2deg.wset.avg456"
case2=case1
dir1="h0zm.seasavg.1960-1979"
dir2="h0zm.seasavg.1986-2005"

season="DJF"
;season="JJA"

;dir1="h0zm.seasavg.1960-1979"
;dir1="h0zm.seasavg.2005-2024"

infile1=case1+".cam2.h0zm."+season+".nc"
infile2=case2+".cam2.h0zm."+season+".nc"

file1=addfile(rootdir+"/"+case1+"/"+dir1+"/"+infile1,"r")
file2=addfile(rootdir+"/"+case2+"/"+dir2+"/"+infile2,"r")

U1=file1->U
U2=file2->U
lev=U1&lev
lev@units="hPa"
lev@long_name="Pressure"
U1&lev=lev
U2&lev=lev
Udiff=U2 ; get metadata
Udiff=Udiff-U1
;Udiff!0=U1!0
;Udiff!1=U1!1

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
  
  res2 = True
  res2@gsnDraw               = False
  res2@gsnFrame              = False
  res2@cnLinesOn             = True
  res2@cnFillOn              = False

  imagefile="$IMAGE_PATH/"+case2+"_Ujet_"+season
  wks = gsn_open_wks("ps",imagefile)
  gsn_define_colormap(wks,"blwhre")                ; choose colormap
  
  minP = 10
  maxP = 1013.
  
  plot1 = gsn_csm_pres_hgt(wks, Udiff({minP:maxP}, {-90.0:90.0}), res)
  plot2 = gsn_csm_pres_hgt(wks, U1   ({minP:maxP}, {-90.0:90.0}), res2)
  overlay(plot1, plot2)
  draw(plot1)
  frame(wks)
  delete(wks)
  
  muttcmd="mutt mikejmills@mac.com -a "+imagefile+".ps -s "+imagefile+" < /dev/null"
  print(muttcmd)
  mutter=systemfunc(muttcmd)
  print(mutter)
