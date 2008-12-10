; ***********************************************
; sizeDistEnsemble.ncl
; ***********************************************
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "tests/core.ncl"
load "tests/getHistoryPath.ncl"
;************************************************

undef("sizeDistEnsemble")
procedure sizeDistEnsemble(caseid:string, targetLat:float, targetDate:string, \
                           targetLev:float, file1:string, file2:string)

begin
  print("")
  print(caseid)
  historyPath=getHistoryPath(caseid)
  
  ifile=historyPath+"daily/"+file1
 ; Open the sulfate file
  print("")
  print("  Opening sulfate file ...")
  print("    " + ifile)
  history = addfile(ifile, "r")
  
 ; Open the dust file
  print("")
  print("  Opening dust file ...")
  print("    " + file2)
  dustfile = addfile(file2, "r")

  time  = history->time
  date  = history->date
  
 ; vars(time, lev, lat, lon, bins)
  print("  reading DUSTNDBIN")
  dustndbin  = history->DUSTNDBIN
  print("  reading DUSTRAD")
  dustrad  = history->DUSTRAD

 ; vars(time, lev, lat, lon, bins)
  print("  reading MIXSULNDBIN")
  mixsulndbin  = history->MIXSULNDBIN
  print("  reading SULFRAD")
  sulfrad  = history->SULFRAD
  print("  reading PURSULNDBIN")
  pursulndbin  = history->PURSULNDBIN
  print("  reading PURERAD")
  purerad  = history->PURERAD

 ; vars(time, lev, lat, lon)
  print("  reading Z3")
  z3  = history->Z3

  lat = sulfrad&lat
  nLats = dimsizes(lat)

;    nBins=dimsizes(dustrad)

 ; ratio of mass in bin i to bin i+1:
  rmrat = 2.0
  dlog10r = (1./3.)*log10(rmrat)

  iTime=ind(date.eq.targetDate)

 ;************************************************
 ; to plot multiple lines, you must put them into
 ; a mulidimensional array
 ;************************************************

  data   = mixsulndbin(iTime,{targetLev},{targetLat},:,:)
  radius = sulfrad(iTime,{targetLev},{targetLat},:,:)
  
  data2   = pursulndbin(iTime,{targetLev},{targetLat},:,:)
  radius2 = purerad(iTime,{targetLev},{targetLat},:,:)
  
  data3   = dustndbin(iTime,{targetLev},{targetLat},:,:)
  radius3 = dustrad(:)

 ; convert to dN/dlog10r
  data = data / dlog10r
  data2 = data2 / dlog10r
  data3 = data3 / dlog10r
  
  print("mixsulf:"+min(data)+" | "+max(data))
  print("pursulf:"+min(data2)+" | "+max(data2))
  print("dust:"+min(data3)+" | "+max(data3))

  data@long_name="dN/dlogr (# cm~S~-3~N~)"
  data@units="/cm3"
  data@_FillValue=0.0
  
  data2@long_name="dN/dlogr (# cm~S~-3~N~)"
  data2@units="/cm3"
  data2@_FillValue=0.0
  
  data3@long_name="dN/dlogr (# cm~S~-3~N~)"
  data3@units="/cm3"
  data3@_FillValue=0.0

  radius=radius*1e7  ; convert cm -> nm
  radius@long_name="sulfate radius (nm)"
  radius@units="nm"
  
  radius2=radius2*1e7  ; convert cm -> nm
  radius2@long_name="sulfate radius (nm)"
  radius2@units="nm"
  
  radius3=radius3*1e7  ; convert cm -> nm
  radius3@long_name="dust radius (nm)"
  radius3@units="nm"
  
  print("mixsulf:"+min(radius) +" | "+max(radius))
  print("pursulf:"+min(radius2)+" | "+max(radius2))
  print("dust:   "+min(radius3)+" | "+max(radius3))

 ;************************************************
 ; plotting parameters
 ;************************************************
 ; get average geometric altitude to 0.1 km
  zint=floattoint(avg(z3(iTime,{targetLev},{targetLat},:))/100.0+0.5)
  z=int2flt(zint)
;  print(z)
  z=z/10.
;  print(z)

  wks   = gsn_open_wks ("ps","sizeDist."+caseid+"_"+targetDate+"_"+targetLat+"_"+z) ; open workstation

  res                   = True                      ; plot mods desired
  res@gsnDraw           = False
  res@gsnFrame          = False
  res@tiMainString      = ""+targetDate+", "+targetLat+"~S~o~N~N, "+z+" km"       ; add title

 ; note a similiar resource is xyLineThicknessF, which will effect all
 ; the lines in the array.
;    res@xyLineThicknesses = (/6.0,6.0/)               ; make 2nd lines thicker

  res@xyLineColor       = "blue"          ; change line color
  res@xyMonoLineColor   = True            ; same color for all lines
  res@xyDashPattern     = 0                         ; all solid lines
  res@trXLog = True 
  res@trYLog = True 

  res@trYMinF = 1e-5
  res@trYMaxF = 1e5

  res@trXMinF = 0.2  ; 0.2 to 100 nm range x-axis
  res@trXMaxF = 100

;    res@tmEqualizeXYSizes = True
;    res@tmXTLabelFont = 25 ; change the font to Times Roman
  res@tmXTLabelFontHeightF = 0.04 ; Change the font size (default 0.02)

;    res@tiMainFont = 25 ; change the font to Times Roman
;    res@tiMainFontHeightF = 0.04 ; change the font size
;    res@tiUseMainAttributes = True ; use main tile font attributes for X and Y titles

  res@vpXF=0.3  ; move left edge of viewport overs
    
  res2 = res
  res2@xyLineColor = "red"
    
  res3 = res
  res3@xyLineColor = "black"
  
;  nPlots = 2
;  plot = new((/ nPlots /), graphic)

  plot2 = gsn_csm_xy (wks,radius2,data2,res2) ; draw second plot
  plot3 = gsn_csm_xy (wks,radius,data,res)    ; draw third plot
  
  plot  = gsn_csm_xy (wks,radius3,data3,res3) ; create plot
;  overlay(plot, plot2)
  
  resP                      = True
  resP@gsnMaximize          = True
;  resP@gsnPanelLabelBar     = True
;  resP@lbTitleOn            =  True                ; turn on title
;  resP@lbTitlePosition      = "Right"              ; title position
;  resP@lbTitleFontHeightF   = .015                ; make title smaller
;  resP@lbTitleDirection     = "Across"             ; title direction
  resP@gsnPanelXWhiteSpacePercent = 5
  resP@tiMainString          = caseid
  overlay(plot, plot2)
  overlay(plot, plot3)
  
  gsn_panel(wks, plot, (/ 1, 1 /), resP)

end ; procedure sizeDistEnsemble()

targetLat = 78
;targetLat = 82
;targetLev=0.002404575    ; ~86 km at 78N (?)
targetLev=0.00230591145   ; 86.7 km at 78N
;targetLev=0.0016106433  ; 88.1 km at 78N
;targetLev=0.0011211287  ; 89.5 km at 78N
;targetLev=0.00584271015  ; 83 km at 78N
;targetLev=0.00395005   ; 84.3 km at 82N

caseid = "125smin"

;targetDate = "19970701"
;file1 = "all.125smin.cam2.h3.1997-06-30-00000.nc"
;targetDate = "19970709"
;;file1 = "all.125smin.cam2.h3.1997-07-05-00000.nc"
file2 = "/ptmp/bardeen/waccm/carma/run/125_lt_gw2_28/125_lt_gw2_28.cam2.h0.2028-06-15-00000.nc"

;targetDate = "19970715"
;file1 = "all.125smin.cam2.h3.1997-07-15-00000.nc"

;targetDate = "19970615"
;file1 = "all.125smin.cam2.h3.1997-06-15-00000.nc"

targetDate = "19970601"
file1 = "all.125smin.cam2.h3.1997-05-31-00000.nc"
file2 = "/ptmp/bardeen/waccm/carma/run/125_lt_gw2_28/125_lt_gw2_28.cam2.h0.2028-05-16-00000.nc"

;targetDate = "19970801"
;file1 = "all.125smin.cam2.h3.1997-07-30-00000.nc"

sizeDistEnsemble(caseid, targetLat, targetDate, targetLev, file1, file2)
