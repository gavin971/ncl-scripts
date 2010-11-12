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
caraer_rat=caraer_rat/(caraer2+1e-20)
caraer_rat&lev=lev
caraer_rat@units="Experiment/Control"

  wks=gsn_open_wks("pdf","$IMAGE_PATH/bc/CarbonAerRatio")
  gsn_define_colormap(wks,"gui_default")                ; choose colormap
  
  res                       = True
  res@cnFillOn              = True                  ; turn on color fill
  res@gsnSpreadColors       = True                  ; use full range of colors
  res@lbLabelAngleF = 60
  
  plot=gsn_csm_pres_hgt(wks,caraer_rat({0.05:500},:),res)
