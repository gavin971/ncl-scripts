load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

testcase="vgsoot08"
cntlcase="vgcontrol02"

rootpath="$MODEL_PATH/waccm/run/carma-soot/"
climopath="climo/2035-2064/"

testpath=rootpath+testcase+"/"+climopath
testfile=testcase+"_JJA_climo.nc"

cntlpath=rootpath+cntlcase+"/"+climopath
cntlfile=cntlcase+"_JJA_climo.nc"

testf=addfile(testpath+testfile,"r")
cntlf=addfile(cntlpath+cntlfile,"r")

emisfile="$MODEL_PATH/waccm/carma-soot/emissions/soot_virgin_galactic_base.nc"
emisf=addfile(emisfile,"r")

soot_emis=emisf->SOOT
lev_emis=emisf->lev
lev_emis=lev_emis/1e3  ; Convert m -> km

Tsoot=testf->T
Tcntl=cntlf->T
Zsoot=testf->Z3
Zcntl=cntlf->Z3

soot_mmr=testf->SOOT01

lat=testf->lat
lon=testf->lon
lev=testf->lev
time=testf->time

nlev=dimsizes(lev)

DeltaT=new((/2,nlev/),"float")
DeltaT(0,:)=dim_avg_Wrap(Tsoot(0,:,{45},:))-dim_avg_Wrap(Tcntl(0,:,{45},:))
DeltaT(1,:)=dim_avg_Wrap(Tsoot(0,:,{-45},:))-dim_avg_Wrap(Tcntl(0,:,{-45},:))

soot_mass=new((/2,nlev/),"float")
soot_mass(0,:)=dim_avg_Wrap(soot_mmr(0,:,{45},:))
soot_mass(1,:)=dim_avg_Wrap(soot_mmr(0,:,{-45},:))
soot_mass=1e12*soot_mass  ; Convert kg soot/kg air -> ng soot/kg air

z3=new((/2,nlev/),"float")
z3(0,:)=(dim_avg_Wrap(Zsoot(0,:,{45},:))+dim_avg_Wrap(Zcntl(0,:,{45},:)))/2.0
z3(1,:)=(dim_avg_Wrap(Zsoot(0,:,{-45},:))+dim_avg_Wrap(Zcntl(0,:,{-45},:)))/2.0
z3=z3/1e3  ; Convert m -> km

wks   = gsn_open_wks ("pdf", "$IMG_PATH/vgProfilePlot")
;wks   = gsn_open_wks ("pdf", "/Volumes/Data/Images/vgProfilePlot")


  plot = new(4, graphic)

  res = True
  res@gsnDraw           = False
  res@gsnFrame          = False
  
  res@vpXF             = 0.1                  ; position and sizes
  res@vpYF             = 0.94                  ; for XY plot
  res@vpWidthF         = 0.2
  res@vpHeightF        = 0.25
  
  res@trYMinF = 0.0
  res@trYMaxF = 50.0
  res@trXMaxF = 105.0
  res@tiYAxisString = "Altitude (km)"
  res@tiXAxisString = "Soot emission rate (kg m~S~-1~N~ s~S~-1~N~)"

  plot(0) = gsn_csm_xy (wks,soot_emis,lev_emis,res)
  
  res@tiYAxisString = " "
  res@tmYLLabelsOn = False
  res@vpXF             = 0.30                  ; position and sizes
  res@trXMaxF = 8.5
  res@tiXAxisString = "Soot mass (ng/kg air)"
  
  plot(1) = gsn_csm_xy (wks,soot_mass,z3,res)

  res@vpXF             = 0.50                  ; position and sizes
  res@trXMinF = -0.5
  res@trXMaxF = 0.5
  res@tiXAxisString = "T~B~test~N~-T~B~control~N~ (K)"
  
  plot(2) = gsn_csm_xy (wks,DeltaT,z3,res)
  
  res@vpXF             = 0.70                  ; position and sizes
  res@tmYROn = True
  res@tmYRLabelsOn = True
  plot(3) = gsn_csm_xy (wks,DeltaT,z3,res)

;
; Maximize plots in frame.
;
 psres = True   
; psres@gsnPaperOrientation = "landscape"                                                   
 maximize_output(wks,psres)  ; calls draw and frame for you 

;
; If you try to use gsn_panel to panel different sized plots, you
; will get unexpected results. You can turn on gsnPanelDebug to
; see what viewport values are used, and then use these values
; to set the viewport resources for each individual plot.
;  resP                      = True
;  resP@gsnMaximize          = True
;  resP@gsnPanelXWhiteSpacePercent = 0
;  resP@gsnPanelYWhiteSpacePercent = 0
;  gsn_panel(wks, plot, (/ 1, 2 /), resP)

delete(wks)
