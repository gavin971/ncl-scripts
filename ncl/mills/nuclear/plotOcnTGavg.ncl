load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

host   = systemfunc("hostname")
ostype = systemfunc("echo $OSTYPE")
MODEL_PATH=getenv("MODEL_PATH")
IMAGE_PATH=getenv("IMG_PATH")

caseid = (/ \
   "b40.rcp4_5.2deg.wcm.carma.bc0tg.006", \
   "b40.rcp4_5.2deg.wcm.carma.bc0tg.007", \
   "b40.rcp4_5.2deg.wcm.carma.bc0tg.008", \
   "b40.rcp4_5.2deg.wcm.carma.bc5tgIP.003", \
   "b40.rcp4_5.2deg.wcm.carma.bc5tgIP.006", \
   "b40.rcp4_5.2deg.wcm.carma.bc5tgIP.007" /)
   
rootpath = "/glade/p/acd/mmills/run/"
   
nCases = dimsizes(caseid)

do i = 0,nCases-1
  filename=rootpath+caseid(i)+"/ocn/proc/"+caseid(i)+".pop.h.TEMP.Gavg.2013-2038.nc"
  print(i+" "+filename)
  fileid=addfile(filename,"r")
  T=fileid->TEMP
  if (i .eq. 0) then
    ds=dimsizes(T)
    nTimes=ds(0)
    nLev=ds(1)
    temp = new((/nCases,nTimes,nLev/),"double")
    temp!1=T!0
    temp!2=T!1
    time=T&time
    z_t=T&z_t
    temp&time=time
    temp&z_t=z_t
  end if
  temp(i,:,:)=T
  delete(fileid)
  delete(T)
  delete(filename)
end do   

; z_t = 500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500, 9500, 10500, 
;    11500, 12500, 13500, 14500, 15500, 16509.84, 17547.9, 18629.13, 19766.03, 
;    20971.14, 22257.83, 23640.88, 25137.02, 26765.42, 28548.37, 30511.92, 
;    32686.8, 35109.35, 37822.76, 40878.46, 44337.77, 48273.67, 52772.8, 
;    57937.29, 63886.26, 70756.33, 78700.25, 87882.52, 98470.59, 110620.4, 
;    124456.7, 140049.7, 157394.6, 176400.3, 196894.4, 218645.7, 241397.2, 
;    264900.1, 288938.5, 313340.5, 337979.3, 362767, 387645.2, 412576.8, 
;    437539.2, 462519, 487508.3, 512502.8, 537500 ;

levels=(/5.,55.,105.,200.,300.,500.,1000./) ; levels in meters
levels@units="meters"
nLevels=dimsizes(levels)
nl3=nLevels*3

deltaT=new((/nl3,nTimes/),"double")
deltaT!1=temp!1
deltaT!0="levels"
deltaT&time=time

do i = 0, nLevels-1
  lev=levels(i)*100.
  cntlavg=dim_avg_n_Wrap(temp(0:2,:,{lev}),0)
  j=nLevels*2+i
  deltaT(j,:)   = dim_avg_n(temp(3:5,:,{lev}),0) - cntlavg
  deltaT(2*i,:)   = dim_min_n(temp(3:5,:,{lev}),0) - cntlavg
  deltaT(2*i+1,:) = dim_max_n(temp(3:5,:,{lev}),0) - cntlavg
end do

pdffile="Fig06_OceanTempVsTime"
wks   = gsn_open_wks ("pdf",IMAGE_PATH+"/nuclear/"+pdffile) ; open workstation

res                   = True                      ; plot mods desired
res@vpWidthF  = 0.7          ; set width of plot
;res@vpHeightF = 1.0          ; set height of plot
res@xyMonoLineColor   = False   ; different colors for different lines
res@xyMonoDashPattern = True 
res@xyLineColors      = (/ \
    "grey","grey","plum","plum","lightblue","lightblue","darkolivegreen1","darkolivegreen1", \
    "gold","gold","pink","pink","darkgoldenrod1","darkgoldenrod1",\
    "black","purple4","blue","darkgreen","orange3","red","sienna4" /)
res@gsnXYFillColors = (/\
    "grey","transparent","plum","transparent","lightblue","transparent","darkolivegreen1","transparent", \
    "gold","transparent","pink","transparent","darkgoldenrod1","transparent",\
    "transparent","transparent","transparent","transparent","transparent","transparent" /)
res@xyLineThicknessF = 3.0
res@xyMonoLineThickness = False
res@xyLineThicknesses = (/0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5, \
    3.0,3.0,3.0,3.0,3.0,3.0,3.0/)
res@trXMaxF=26.0
res@trYMaxF=0.1
res@trYMinF=-0.9
res@tiXAxisString="Time (years since BC injection)"
res@tiYAxisString="Temperature (K)"
res@gsnFrame  = False
res@gsnPaperOrientation    = "landscape"
res@tmXBMode="Manual"
res@tmXBTickSpacingF=2
res@tmXBAutoPrecision=False
res@tmXBPrecision=2

time=(time-time(0))/365

plot = gsn_csm_xy (wks,time,deltaT,res) ; draw plot

txres = True
txres@txFontHeightF = 0.03
do i = 0, 3
  txres@txFontColor=res@xyLineColors(nLevels*2+i)
  gsn_text_ndc(wks,floattoint(levels(i))+" m",0.8,0.25+.065*i,txres)
end do
do i = 4, 6
  txres@txFontColor=res@xyLineColors(nLevels*2+i)
  gsn_text_ndc(wks,floattoint(levels(i))+" m",0.52,0.32+.072*i,txres)
end do

frame(wks)

delete(wks)
