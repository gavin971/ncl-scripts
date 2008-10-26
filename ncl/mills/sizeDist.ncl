; ***********************************************
; sizeDist.ncl
; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;************************************************
begin

 historyRun="geo2"
 histFile="geo2.cam2.h0.1996-07.nc"
 controlRun="nongeo2"
 contFile="nongeo2.cam2.h0.1996-07.nc"
 
 iLev=0.07991075
 iLat=46

 sizehist=addfile(historyRun+"/"+histFile,"r")
 sizecont=addfile(controlRun+"/"+contFile,"r")

 time=sizehist->time
 lev=sizehist->lev
 lat=sizehist->lat
 lon=sizehist->lon

 mixsulndbin1=sizehist->MIXSULNDBIN ; ( time, lev, lat, lon, bins )
 sulfrad1=sizehist->SULFRAD  ; ( time, lev, lat, lon, bins )

 mixsulndbin2=sizehist->MIXSULNDBIN ; ( time, lev, lat, lon, bins )
 sulfrad2=sizehist->SULFRAD  ; ( time, lev, lat, lon, bins )

; convert cm -> microns
 sulfrad1zavg=1e4*dim_avg_Wrap(sulfrad1(time|0,lev|{iLev},lat|{iLat},bins|:,lon|:))
 sulfrad2zavg=1e4*dim_avg_Wrap(sulfrad2(time|0,lev|{iLev},lat|{iLat},bins|:,lon|:))\
 
 mixsulnd1zavg=dim_avg_Wrap(mixsulndbin1(time|0,lev|{iLev},lat|{iLat},bins|:,lon|:))
 mixsulnd2zavg=dim_avg_Wrap(mixsulndbin2(time|0,lev|{iLev},lat|{iLat},bins|:,lon|:))
 
 nBins=38

 mixsuldr1=new(nBins)
 mixsuldr2=new(nBins)
 do i = 1,nBins-2
   mixsuldr1(i)= sqrt(sulfrad1zavg(i+1)* sulfrad1zavg(i)) - sqrt(sulfrad1zavg(i)* sulfrad1zavg(i-1))
   mixsuldr2(i)= sqrt(sulfrad2zavg(i+1)* sulfrad2zavg(i)) - sqrt(sulfrad2zavg(i)* sulfrad2zavg(i-1))
 end do
 
 mixsuldr1(nBins-1) = mixsuldr1(nBins-2) * mixsuldr1(nBins-2)/mixsuldr1(nBins-3)
 mixsuldr2(nBins-1) = mixsuldr2(nBins-2) * mixsuldr2(nBins-2)/mixsuldr2(nBins-3)
 
;************************************************
; to plot multiple lines, you must put them into
; a mulidimensional array
;************************************************
 data      = new((/2,nBins/),float)
 radius    = new((/2,nBins/),float)
 
 radius(0,:) = sulfrad1zavg(:)
 radius(1,:) = sulfrad2zavg(:)
; radius=radius*1000.  ; convert microns -> nm
 radius@long_name="particle radius (nm)"
 radius@units="nm"

; Scale to dN/dlogr
;  dN/d(log10(r) = dN/dr * r * ln(10.)
 data(0,:) = (mixsulnd1zavg(:) / mixsuldr1(:)) * sulfrad1zavg(:) * log(10.)
 data(1,:) = (mixsulnd2zavg(:) / mixsuldr2(:)) * sulfrad2zavg(:) * log(10.)
 
 data@long_name="dN/dlogr (# cm-3)"
 data@units="/cm3"
 data@_FillValue=0.0
 print(max(data))
 print(min(data))
 
 radius=radius*1000.  ; convert microns -> nm
 
;************************************************
; plotting parameters
;************************************************
 wks   = gsn_open_wks ("ps","sizeDist."+historyRun) ; open workstation

 res                   = True                      ; plot mods desired
 res@tiMainString      = ""+iLev+" "+iLat+" "+iLon       ; add title

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
