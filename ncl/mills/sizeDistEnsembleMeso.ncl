; ***********************************************
; sizeDistEnsembleMeso.ncl
; ***********************************************
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "tests/core.ncl"
load "tests/getHistoryPath.ncl"
;************************************************

undef("sizeDistEnsembleMeso")
procedure sizeDistEnsembleMeso(targetLat:float, targetLev:float, caseid:string)

begin

  aFile="analysis/"+caseid+"/sizeDistAnalysis."+caseid+".nc"
  print("  Opening analysis file "+aFile)
  history = addfile(aFile, "r")
  
 ; vars(time, lev, lat, lon, bins)
  lev = history->lev
  lat = history->lat
  date= history->date

  print("  reading mixsulndbin")
  mixsulndbin  = history->mixsulndbin
  print("  reading mixsulrad")
  mixsulrad  = history->mixsulrad
  print("  reading pursulndbin")
  pursulndbin  = history->pursulndbin
  print("  reading pursulrad")
  pursulrad  = history->pursulrad
  
  if (.not.(isfilepresent("images/"+caseid))) then
    a = systemfunc("mkdir -p images/"+caseid)
  end if
  
  wks   = gsn_open_wks ("pdf","images/"+caseid+"/sizeDistEnsembleMeso."+caseid) ; open workstation

  res                   = True                      ; plot mods desired
  res@gsnDraw           = False
  res@gsnFrame          = False
    ; add title

 ; note a similiar resource is xyLineThicknessF, which will effect all
 ; the lines in the array.
;    res@xyLineThicknesses = (/6.0,6.0/)               ; make 2nd lines thicker

  res@xyLineColor       = "blue"          ; change line color
  res@xyMonoLineColor   = True            ; same color for all lines
  res@xyDashPattern     = 0                         ; all solid lines
  res@trXLog = True 
  res@trYLog = True 

  res@trYMinF = 1e-5
  res@trYMaxF = 1e4

  res@trXMinF = 0.2  ; 0.2 to 100 nm range x-axis
  res@trXMaxF = 30

;    res@tmEqualizeXYSizes = True
;    res@tmXTLabelFont = 25 ; change the font to Times Roman
  res@tmXTLabelFontHeightF = 0.04 ; Change the font size (default 0.02)

;    res@tiMainFont = 25 ; change the font to Times Roman
;    res@tiMainFontHeightF = 0.04 ; change the font size
;    res@tiUseMainAttributes = True ; use main tile font attributes for X and Y titles

  res@vpXF=0.3  ; move left edge of viewport overs
    
  res2 = res
  res2@xyLineColor = "red"

  plot = new(4, graphic)
  plot2 = new(4, graphic)
  
  do i = 0,3
    print(date(i))

    res@tiMainString      = ""+date(i)+" "+lat({targetLat})+"~S~o~N~N "+lev({targetLev})+" hPa"
   ; ratio of mass in bin i to bin i+1:
    rmrat = 2.0
    dlog10r = (1./3.)*log10(rmrat)

   ;************************************************
   ; to plot multiple lines, you must put them into
   ; a mulidimensional array
   ;************************************************

    data   = mixsulndbin(i,{targetLev},{targetLat},:,:)
    radius = mixsulrad(i,{targetLev},{targetLat},:,:)

    data2   = pursulndbin(i,{targetLev},{targetLat},:,:)
    radius2 = pursulrad(i,{targetLev},{targetLat},:,:)

   ; convert to dN/dlog10r
    data = data / dlog10r
    data2 = data2 / dlog10r

    print("mixsulf:"+min(data)+" | "+max(data))
    print("pursulf:"+min(data2)+" | "+max(data2))

    data@long_name="dN/dlogr (# cm~S~-3~N~)"
    data@units="/cm3"
    data@_FillValue=0.0

    data2@long_name="dN/dlogr (# cm~S~-3~N~)"
    data2@units="/cm3"
    data2@_FillValue=0.0

    radius=radius*1e7  ; convert cm -> nm
    radius@long_name="sulfate radius (nm)"
    radius@units="nm"

    radius2=radius2*1e7  ; convert cm -> nm
    radius2@long_name="sulfate radius (nm)"
    radius2@units="nm"

    print("mixsulf:"+min(radius) +" | "+max(radius))
    print("pursulf:"+min(radius2)+" | "+max(radius2))

   ;************************************************
   ; plotting parameters
   ;************************************************
   ; get average geometric altitude to 0.1 km


  ;  nPlots = 2
  ;  plot = new((/ nPlots /), graphic)

    plot2(i) = gsn_csm_xy (wks,radius2,data2,res2) ; draw second plot
    plot(i) = gsn_csm_xy (wks,radius,data,res)    ; draw first plot

    resP                      = True
    resP@gsnMaximize          = True
    resP@gsnPanelXWhiteSpacePercent = 5
  ;  resP@tiMainString          = caseid
    overlay(plot(i), plot2(i))
  end do
  
  gsn_panel(wks, plot, (/ 2, 2 /), resP)

end ; procedure sizeDistEnsembleMeso()

targetLat = 78
;targetLat = 82
;targetLev=0.002404575    ; ~86 km at 78N (?)
targetLev=0.00230591145   ; 86.7 km at 78N
;targetLev=0.0016106433  ; 88.1 km at 78N
;targetLev=0.0011211287  ; 89.5 km at 78N
;targetLev=0.00584271015  ; 83 km at 78N
;targetLev=0.00395005   ; 84.3 km at 82N

caseid="125b36dq"

sizeDistEnsembleMeso(targetLat, targetLev, caseid)
