load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("plotvgDiffAnn")
procedure plotvgDiffAnn(imgPath:string, specname:string,\
   species_x[*]:float, species_c[*]:float, \
   speciesSD_x[*]:float, speciesSD_c[*]:float, \
   species2d_x[*][*]:float, species2d_c[*][*]:float, \
   speciesSD2d_x[*][*]:float, speciesSD2d_c[*][*]:float)
   
local species_diff,wks,res,lev
begin 

  lat=species_x&lat
  nLat=(/dimsizes(lat)/)
  
  lev = species2d_x&lev
  
  species_diff=species_x-species_c
  species_diff!0 = "lat"
  species_diff&lat = lat

  species2d_diff=species2d_x-species2d_c
  species2d_diff!0 = "lev"
  species2d_diff&lev = lev
  species2d_diff!1 = "lat"
  species2d_diff&lat = lat

  species_dpct=100.*species_diff/species_c
  species_dpct!0 = "lat"
  species_dpct&lat = lat
  
  species2d_dpct=100.*species2d_diff/species2d_c
  species2d_dpct!0 = "lev"
  species2d_dpct&lev = lev
  species2d_dpct!1 = "lat"
  species2d_dpct&lat = lat

  
  iflag= False               ; population variance similar
  prob = ttest(species_x, speciesSD_x^2.0, 30, species_c, speciesSD_c^2.0, 30, iflag, False)
  prob!0 = "lat"
  prob&lat = lat
  
  prob2d = ttest(species2d_x, speciesSD2d_x^2.0, 30, species2d_c, speciesSD2d_c^2.0, 30, iflag, False)
  prob2d!0 = "lev"
  prob2d&lev = lev
  prob2d!1 = "lat"
  prob2d&lat = lat
  
  sig_level = 0.025  ; significance level (95%)
  
  flag=(prob.lt.sig_level)
  
  fmt = "%5.2f"
  emt = "%10.5e"
  gmt = "%4.0f"
  
  do i = 0,nLat-1
    print(sprintf(gmt,lat(i))+"  "+sprintf(fmt, species_x(i))+"  "+sprintf(fmt, species_c(i))+ \
       "  "+sprintf(fmt, species_diff(i))+"  "+sprintf(fmt, species_dpct(i))+"  "+\
       sprintf(emt, prob(i))+" "+flag(i))
  end do
  
  
  do i = 0,nLat-1
    print(sprintf(gmt,lat(i))+"  "+sprintf(emt, species2d_x({100},i))+"  "+sprintf(emt, species2d_c({100},i))+ \
       "  "+sprintf(emt, speciesSD2d_x({100},i))+"  "+sprintf(emt, speciesSD2d_c({100},i))+"  "+\
       sprintf(emt, prob2d({100},i)))
  end do
  
       
  if (.not.(isfilepresent(imgPath+"/vg"))) then
    a = systemfunc("mkdir -p "+imgPath+"/vg")
  end if  
  wks = gsn_open_wks("pdf", imgPath+"/vg/vg"+specname+"ChangeAnn")
    
  resP                      = True
  resP@gsnMaximize          = True
  resP@gsnPanelXWhiteSpacePercent = 5
  resP@gsnPaperOrientation = "landscape"
  
  res                     = True          ; plot mods desired
  res@gsnPaperOrientation = "landscape"

  res@gsnMaximize         = True          ; expand plot size to maximum
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.7
;  res@tiXAxisString       = "Latitude"

;  res@tmXBLabelFontHeightF =0.025
;  res@tmYLLabelFontHeightF =0.03
  
;===================================================================  
;  res@gsnRightString      = "O~B~3~N~ Column(DU)"
;  res@gsnLeftString       = "Control run"
  
  data=new((/2,nLat/),float)
  data(0,:)=species_c
  data(1,:)=species_x
  gsn_define_colormap(wks,"rainbow")
  plot  = gsn_csm_xy (wks,lat,data,res) ; create plot
  
;  res@gsnDraw             = False       ; Do not draw plot
;  res@gsnFrame            = False       ; Do not advance frame
  data2=new((/5,nLat/),float)
  data2(0,:)=species_diff
  data2(1,:)=species2d_diff({100},:)
  data2(2,:)=species2d_diff({60},:)
  data2(3,:)=species2d_diff({20},:)
  data2(4,:)=0.0
  markers2=new((/4,nLat/),float)
  markers2(0,:)=where(prob.lt.sig_level,data2(0,:),-999)
  markers2(1,:)=where(prob2d({100},:).lt.sig_level,data2(1,:),-999)
;  print("markers2(1):"+markers2(1,:)+" "+prob2d({100},:))
  markers2(2,:)=where(prob2d({60},:).lt.sig_level,data2(2,:),-999)
  markers2(3,:)=where(prob2d({20},:).lt.sig_level,data2(3,:),-999)
  res@gsnDraw             = False       ; Do not draw plot
  res@gsnFrame            = False       ; Do not advance frame
  res@xyMarkLineMode   = "Lines"                ; choose which have markers
  res@xyMonoDashPattern = True
  res@xyLineThicknessF = 1.0
  res2=res
  res2@xyLineThicknessF = 10.0
  res2@xyLineColors = (/ "grey50" ,"dodgerblue4", "firebrick4", "darkolivegreen" /)
  res@xyLineColors = (/ "black", "dodgerblue1", "firebrick1", "darkgreen", "black" /)
  plot  = gsn_csm_xy (wks,lat,markers2,res2) ; create plot
  plot2  = gsn_csm_xy (wks,lat,data2,res) ; create plot
  overlay (plot, plot2)
  draw(plot)
  frame(wks)
  
  res@gsnDraw             = False       ; Do not draw plot
  res@gsnFrame            = False       ; Do not advance frame
  res@xyMarkLineMode   = "Lines" 
;  res@xyMarkLineMode   = "Markers" 
  data3=new((/5,nLat/),float)
  data3(0,:)=species_dpct
  data3(1,:)=species2d_dpct({100},:)
  data3(2,:)=species2d_dpct({60},:)
  data3(3,:)=species2d_dpct({20},:)
  data3(4,:)=0.0
  markers3=new((/4,nLat/),float)
  markers3(0,:)=where(prob.lt.sig_level,data3(0,:),-999)
  markers3(1,:)=where(prob2d({100},:).lt.sig_level,data3(1,:),-999)
  markers3(2,:)=where(prob2d({60},:).lt.sig_level,data3(2,:),-999)
  markers3(3,:)=where(prob2d({20},:).lt.sig_level,data3(3,:),-999)
;  res@xyMarkers   = markers2             ; choose which have markers
  plot  = gsn_csm_xy (wks,lat,markers3,res2) ; create plot
  plot2  = gsn_csm_xy (wks,lat,data3,res) ; create plot
  overlay (plot, plot2)
  draw(plot)
  frame(wks)
  
end


undef("plotvgChangeAnn")
procedure plotvgChangeAnn(rootPath:string,imgPath:string,historyRun:string, controlRun:string)
local T_x,T_c

begin
  print((/historyRun/))
  sootAnnMeanFile = rootPath+"/"+historyRun+"/h0zm/"+historyRun+".cam2.h0zm.annual.nc"
  print((/sootAnnMeanFile/))
  sootAnnSDFile = rootPath+"/"+historyRun+"/h0zm/"+historyRun+".cam2.h0zm.annual_sd.nc"
  print((/sootAnnSDFile/))
  
  print((/controlRun/))
  controlAnnMeanFile = rootPath+"/"+controlRun+"/h0zm/"+controlRun+".cam2.h0zm.annual.nc"
  print((/controlAnnMeanFile/))
  controlAnnSDFile = rootPath+"/"+controlRun+"/h0zm/"+controlRun+".cam2.h0zm.annual_sd.nc"
  print((/controlAnnSDFile/))
  
  a=addfile(sootAnnMeanFile,"r")
    O3_col_x=a->O3_col
    TS_x=a->TS
    O3_x=a->O3
    T_x=a->T
  b=addfile(sootAnnSDFile,"r")
    O3_colSD_x=b->O3_col(0,:)
    TSSD_x=b->TS(0,:)
    O3SD_x=b->O3(0,:,:)
    TSD_x=b->T(0,:,:)
  c=addfile(controlAnnMeanFile,"r")
    O3_col_c=c->O3_col
    TS_c=c->TS
    O3_c=c->O3
    T_c=c->T
  d=addfile(controlAnnSDFile,"r")
    O3_colSD_c=d->O3_col(0,:)
    TSSD_c=d->TS(0,:)
    O3SD_c=d->O3(0,:,:)
    TSD_c=d->T(0,:,:)
 
  print("plotting ozone column...")
  plotvgDiffAnn(imgPath, "O3", O3_col_x, O3_col_c, O3_colSD_x, O3_colSD_c,\
          O3_x, O3_c, O3SD_x, O3SD_c)
  plotvgDiffAnn(imgPath, "Temp", TS_x, TS_c, TSSD_x, TSSD_c,\
          T_x, T_c, TSD_x, TSD_c)
  
end

historyRun="vgsoot08"
controlRun="vgcontrol02"
;rootPath="/data5/mmills"
;rootPath="/Volumes/Data/Models/waccm/carma-soot/run"
rootPath="/Data/Models/waccm/run/carma-soot"
imgPath="/Data/Models/waccm/images"

plotvgChangeAnn(rootPath,imgPath,historyRun,controlRun)
