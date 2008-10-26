; ***********************************************
; sizeDist.ncl
; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;************************************************
begin

 historyRun="JPL06_ncar"
 histFile="all.JPL06_1.cam2.h0.1995-07-25-00000.nc"
 iTime=208
; histFile="all.JPL06_1.cam2.h0.1995-05-31-00000.nc"
; iTime=150
 sizeFile=histFile
 ;history=addfile(historyRun+"/"+histFile,"r")
 sizehist=addfile(historyRun+"/"+sizeFile,"r")

 time=sizehist->time
 lev=sizehist->lev
 lat=sizehist->lat
 lon=sizehist->lon

 dustndbin=sizehist->DUSTNDBIN ; ( time, lev, lat, lon, bins )
 dustrad=sizehist->DUSTRAD ; ( bins )
 dustdr=sizehist->DUSTDR ; ( bins )
 mixsulndbin=sizehist->MIXSULNDBIN ; ( time, lev, lat, lon, bins )
 sulfrad=sizehist->SULFRAD  ; ( time, lev, lat, lon, bins )

; convert cm -> microns
 dustrad=dustrad*1e4
 dustdr=dustdr*1e4
 sulfrad=sulfrad*1e4
 
 nBins=dimsizes(dustrad)
 
 iLev=0.002404575
 iLat=86
 iLon=75
 
; ratio of mass in bin i to bin i+1:
 rmrat = 2.0
 
;
 
 mixsuldr=dustdr
 do i = 1,nBins-2
;   mixsuldr(i)= sqrt(sulfrad({iTime},{iLev},{iLat},{iLon},i+1)* sulfrad({iTime},{iLev},{iLat},{iLon},i)) - sqrt(sulfrad({iTime},{iLev},{iLat},{iLon},i)* sulfrad({iTime},{iLev},{iLat},{iLon},i-1))
   mixsuldr(i) = 
 end do
 
 mixsuldr(nBins-1) = mixsuldr(nBins-2) * mixsuldr(nBins-2)/mixsuldr(nBins-3)
 
; do i = 0,nBins-1
;   print("i:"+i+" dusr:"+dustrad(i)+" dusdr:"+dustdr(i)+" sulr:"+sulfrad({iTime},{iLev},{iLat},{iLon},i)+" ;suldr:"+mixsuldr(i))
; end do

;************************************************
; to plot multiple lines, you must put them into
; a mulidimensional array
;************************************************
 data      = new((/2,nBins/),float)
 radius    = new((/2,nBins/),float)
 
 radius(0,:) = dustrad(:)
 radius(1,:) = sulfrad({iTime},{iLev},{iLat},{iLon},:)
; radius=radius*1000.  ; convert microns -> nm
 radius@long_name="particle radius (nm)"
 radius@units="nm"

; Scale to dN/dlog10r
;  dN/d(log10(r) = dN/dr * r * ln(10.)
 data(0,:) = (dustndbin({iTime},{iLev},{iLat},{iLon},:) / dustdr(:)) * dustrad(:) * log(10.)
 data(1,:) = (mixsulndbin({iTime},{iLev},{iLat},{iLon},:) /mixsuldr(:)) * radius(1,:)  * log(10.)
 
 data@long_name="dN/dlogr (# cm-3)"
 data@units="/cm3"
 data@_FillValue=0.0
 print(max(data))
 print(min(data))
 
 radius=radius*1000.  ; convert microns -> nm
 
;************************************************
; plotting parameters
;************************************************
 wks   = gsn_open_wks ("ps","sizeDist."+historyRun+"."+iTime) ; open workstation

 res                   = True                      ; plot mods desired
 res@tiMainString      = ""+iTime+" "+iLev+" "+iLat+" "+iLon       ; add title

; note a similiar resource is xyLineThicknessF, which will effect all
; the lines in the array.
 res@xyLineThicknesses = (/6.0,6.0/)               ; make 2nd lines thicker

 res@xyLineColors      = (/"blue","red"/)          ; change line color
 res@xyDashPattern     = 0                         ; all solid lines
 res@trXLog = True 
 res@trYLog = True 
 
 res@trYMinF = 0.01
 res@trYMaxF = 1e5
 
 res@trXMinF = 0.1  ; 0.1 to 100 nm range x-axis
 res@trXMaxF = 100
 
 res@tmEqualizeXYSizes = True
 res@tmXTLabelFont = 25 ; change the font to Times Roman
 res@tmXTLabelFontHeightF = 0.04 ; Change the font size (default 0.02)
 
 res@tiMainFont = 25 ; change the font to Times Roman
 res@tiMainFontHeightF = 0.04 ; change the font size
 res@tiUseMainAttributes = True ; use main tile font attributes for X and Y titles
 
 res@vpXF=0.3  ; move left edge of viewport overs
 
 plot  = gsn_csm_xy (wks,radius,data,res) ; create plot
 
end
