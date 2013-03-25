load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

host   = systemfunc("hostname")
ostype = systemfunc("echo $OSTYPE")
MODEL_PATH=getenv("MODEL_PATH")
IMAGE_PATH=getenv("IMG_PATH")

if (host.eq."modeling1") then
  path="/data4/mmills/"
else if (isStrSubset(host,"caldera").or.isStrSubset(host,"geyser").or.isStrSubset(host,"yslogin")) then
  path=MODEL_PATH+"/"
else
  path=MODEL_PATH+"/cesm/case/"  
end if
end if
;path="~/csm/"
print("path="+path)

root="b40.rcp4_5.2deg.wcm.carma.bc0tg"
sub="h0gm"

      
avgfile=path+root+".avg/"+sub+"/"+root+".cam.h0.avg.nc"
print("adding avg file: "+avgfile)
avgfileh=addfile(avgfile,"r")

gw   = avgfileh->gw
lat  = avgfileh->lat
time = avgfileh->time

time   = time/365.0
nLats  = dimsizes(lat)
nTimes = dimsizes(time)
print("nTimes="+nTimes)

; Net shortwave flux at the surface
print("reading FSNS")
FSNS = avgfileh->FSNS(0:nTimes-1)

; Surface temperature
print("reading TS")
TS = avgfileh->TS(0:nTimes-1)

; Total precipitation
print("reading PRECT")
PRECT = avgfileh->PRECT(0:nTimes-1)

; convert m/s to mm/day
PRECT = PRECT * 24. * 3600. * 1000.
 
  pdffile="TS_Control."+root
  wks   = gsn_open_wks ("pdf",IMAGE_PATH+"/"+pdffile) ; open workstation

  nl = integertochar(10)  ; newline character

  ; Create a panel plot
  resP = True
  nPanel=2
  plot = new(nPanel, graphic)
    
  res                   = True                      ; plot mods desired
  res@tmYMajorGrid = True
  res@tmYMajorGridLineDashPattern = 2
  res@vpWidthF            = 1.3           ; change aspect ratio of plot
  res@vpHeightF           = 0.3
  res@gsnDraw = False
  res@gsnFrame     = False                        ; don't advance frame yet
  res@xyLineThicknessF = 5.0
;  res@vpWidthF=0.85

  res@xyLineColor       = "midnightblue"          ; change line color
  res@xyMonoLineColor   = True            ; same color for all lines

  res@tmXBMaxTicks     = 15
  res@tmXBMinorPerMajor= 3
  res@tmXTMinorPerMajor= 3
  res@tmXBAutoPrecision= False
  res@tmXBPrecision    = 2
  
  res@trXMinF = 0.0
  res@trXMaxF = max(time)
  res@tiXAxisString = "Year"
;  res@tmXTLabelsOn = True
;  res@tmYRLabelsOn = True

;  res@trYMinF = -0.01
;  res@trYMaxF = 5.0
;  res@tmYLPrecision = 1
;  res@tmYLAutoPrecision = False
  
  res@tiYAxisString = "SW Flux (W m~S~-2~N~)"
;  res@tmYLPrecision=2
;  res@tmYLMode = "Manual"
;  res@tmYLTickStartF = -20.0
;  res@tmYLTickSpacingF = 5.0
;  res@tmYLMinorPerMajor = 4
    
  plot(0) = gsn_csm_xy (wks,time,FSNS,res) ; draw plot
  
  res@tiYAxisString = "Temperature (K)"
;  res@trYMaxF = 0.1
;  res@trYMinF = -1.9
;  res@tmYLMode = "Manual"
;  res@tmYLTickStartF = -2.0
;  res@tmYLTickSpacingF = 0.5

  baseplot = gsn_csm_xy (wks,time,TS,res) ; draw plot  
  
;  res@trYMaxF = 0.05
;  res@trYMinF = -0.32
;  res@tmYLPrecision = 1
  res@tiYAxisString = "Precipitation (mm d~S~-1~N~)"
;  res@tmYLMode = "Automatic"

  plot(1) = gsn_csm_xy (wks,time,PRECT,res) ; draw plot
  
  resP@txString = "Years after conflict"
  resP@gsnAttachPlotsXAxis =True
  resP@gsnPaperOrientation = "landscape"
  resP@gsnMaximize         = True          ; expand plot size to maximum
  resP@tmYMajorGrid = True
  resP@tiYAxisFontHeightF=0.015
  newplot = gsn_attach_plots(baseplot,plot,res,resP)
  draw(baseplot)
      
  delete(wks)
  
  if (host.eq."modeling1") then
    muttcmd="mutt mikejmills@mac.com -a "+IMAGE_PATH+"/"+pdffile+".pdf -s "+pdffile+" < /dev/null"
    print(muttcmd)
    mutter=systemfunc(muttcmd)
    print(mutter)
  end if

  if (ostype.eq."darwin") then
    opencmd="open "+IMAGE_PATH+"/"+pdffile+".pdf"
    print(opencmd)
    open=systemfunc(opencmd)
    print(open)
  end if
  
