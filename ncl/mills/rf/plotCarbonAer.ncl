load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"



basedir="/Volumes/Data/Models/waccm/run/carma-soot/"

caseid1="bcrckt03"
rootdir1=basedir+caseid1+"/h0zm.annual/"
file1=caseid1+".cam2.h0zm.aer.annual.nc"

caseid2="bccntl02"
rootdir2=basedir+caseid2+"/h0zm.annual/"
file2=caseid2+".cam2.h0zm.aer.annual.nc"

bcrckt=addfile(rootdir1+file1,"r")
bccntl=addfile(rootdir2+file2,"r")

caraer1=bcrckt->CARAER
caraer2=bccntl->CARAER
lev=bcrckt->lev
lev@units="hPa"
lev@long_name="Pressure (hPa)"

caraer_rat=caraer1  ; get attributes
caraer_rat=100.0*(caraer_rat-caraer2)/(caraer2+1e-20)
caraer_rat&lev=lev

  wks=gsn_open_wks("pdf","$IMAGE_PATH/bc/CarbonAerRatio")
;  gsn_define_colormap(wks,"gui_default")                ; choose colormap
;  gsn_define_colormap(wks,"rainbow")
;  gsn_define_colormap(wks,"BlRe")
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  
  res                       = True
  res@gsnPaperOrientation = "landscape"
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.7
  res@cnFillOn              = True                  ; turn on color fill
  res@gsnSpreadColors       = True                  ; use full range of colors
;  res@lbLabelAngleF = 60
  res@lbTitlePosition       = "Right"              ; title position
  res@lbTitleString = "%"
  res@lbLeftMarginF = 0.00
  res@lbRightMarginF = 0.00
  res@lbAutoManage=False
  res@lbTitleFontHeightF=0.04
  res@lbLabelFontHeightF=0.03
  res@gsnLeftString=" "
  res@gsnRightString=" "
  res@tiYAxisFontHeightF=0.025
  res@tmYLMode="Explicit"
  res@tmYLValues=(/250,100,50,25,10,5,2.5,1,0.5,0.25,0.1,0.05/)
  res@tmYLLabels=(/"250","100","50","25","10","5","2.5","1","0.5","0.25","0.1","0.05"/)
;  res@tmYRMode="Explicit"
;  res@tmYRValues=(/10,20,30,40,50,60,70/)
;  res@tmYRLabels=(/"10","20","30","40","50","60","70"/)
  
  plot=gsn_csm_pres_hgt(wks,caraer_rat({0.025:300},:),res)
