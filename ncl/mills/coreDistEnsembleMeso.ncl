; ***********************************************
; sizeDistEnsembleMeso.ncl
; ***********************************************
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "ncl-scripts/core.ncl"
;************************************************
RHO_DUST = 2.0
PI      = 4.0 * atan(1.0)         ; pi
cpi = 4./3. * PI

undef("coreDistEnsembleMeso")
procedure coreDistEnsembleMeso(targetLat:float, targetLev:float, caseid:string)

begin

  aFile="analysis/"+caseid+"/sizeDistAnalysis."+caseid+".nc"
  print("  Opening analysis file "+aFile)
  history = addfile(aFile, "r")
  
 ; vars(time, lev, lat, lon, bins)
  lev = history->lev
  lat = history->lat
  lon = history->lon
  nLon = dimsizes(lon)
  date= history->date
  print(lat)
  print(lev)
  
  print("  reading mixsulndbin")
  mixsulndbin  = history->mixsulndbin
  
  print("  reading mixsulrad")
  mixsulrad  = history->mixsulrad
  
  print("  reading rmass")
  rmass  = history->rmass
  nBins = dimsizes(rmass)
    
  print("  reading dussulbin")
  dussulbin  = history->dussulbin
  
  if (.not.(isfilepresent("images/"+caseid))) then
    a = systemfunc("mkdir -p images/"+caseid)
  end if
  
  wks   = gsn_open_wks ("pdf","images/"+caseid+"/coreDistEnsembleMeso."+caseid+"."+\
          lat({targetLat})+"N_"+lev({targetLev})+"hPa") ; open workstation

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
  res@trYLog = False 

  res@trYMinF = 0.92
  res@trYMaxF = 1.01

  res@trXMinF = 0.4  ; 0.2 to 100 nm range x-axis
  res@trXMaxF = 20

;    res@tmEqualizeXYSizes = True
;    res@tmXTLabelFont = 25 ; change the font to Times Roman
  res@tmXTLabelFontHeightF = 0.08 ; Change the font size (default 0.02)
  res@tmXBMode = "Explicit"
  res@tmXBValues = (/0.4,0.6,0.8,1.,2.,4.,6.,8.,10.,20./)
  res@tmXBLabels = (/"0.4",".6"," ","1","2","4","6","8","10","20"/)
  res@tiXAxisString = "radius (nm)"
  
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

;    res@tiMainString      = ""+date(i)+" "+lat({targetLat})+"~S~o~N~N "+lev({targetLev})+" hPa"

   ; ratio of mass in bin i to bin i+1:
    rmrat = 2.0
    dlog10r = (1./3.)*log10(rmrat)

   ;************************************************
   ; to plot multiple lines, you must put them into
   ; a mulidimensional array
   ;************************************************

    radius = mixsulrad(i,{targetLev},{targetLat},:,:)
    
    data   = new((/nLon, nBins/), double)
    numcon=mixsulndbin(i,{targetLev},{targetLat},:,:)
    numcon@_FillValue=0.0
    do ibin = 0, nBins - 1
      data(:,ibin) = dussulbin(i,{targetLev},{targetLat},:,ibin)/ \
        (numcon(:,ibin) * rmass(ibin))
      print(ibin+" "+radius(0,ibin)+" "+data(0,ibin))
    end do

    data2   = new((/nLon, nBins/), double)
    do ibin = 0, nBins - 1
      data2(:,ibin) = (dussulbin(i,{targetLev},{targetLat},:,ibin)/ \
        (numcon(:,ibin) * RHO_DUST * cpi))^(1./3.)
      data2(:,ibin) = data2(:,ibin) / mixsulrad(i,{targetLev},{targetLat},:,ibin)
      print(ibin+" "+radius(0,ibin)+" "+data2(0,ibin))
    end do

    print("coremass:"+min(data)+" | "+max(data))
    print("corerad:"+min(data2)+" | "+max(data2))

    data@long_name="core mass ratio"
    data@units=""
    data@_FillValue=0.0

    data2@long_name="core radius ratio"
    data2@units=""
    data2@_FillValue=0.0

    radius=radius*1e7  ; convert cm -> nm
    radius@long_name="sulfate radius (nm)"
    radius@units="nm"

    print("radius:"+min(radius) +" | "+max(radius))

   ;************************************************
   ; plotting parameters
   ;************************************************
   ; get average geometric altitude to 0.1 km


  ;  nPlots = 2
  ;  plot = new((/ nPlots /), graphic)

    plot2(i) = gsn_csm_xy (wks,radius,data2,res2) ; draw second plot
    plot(i) = gsn_csm_xy (wks,radius,data,res)    ; draw first plot

    resP                      = True
    resP@gsnMaximize          = True
    resP@gsnPanelXWhiteSpacePercent = 5
  ;  resP@tiMainString          = caseid
;    overlay(plot(i), plot2(i))
  end do
  
  gsn_panel(wks, plot, (/ 2, 2 /), resP)

end ; procedure sizeDistEnsembleMeso()

;targetLat = 78
targetLat = 86
;targetLat = 82
;targetLev=0.002404575    ; ~86 km at 78N (?)
;targetLev=0.00230591145   ; 86.7 km at 78N
;targetLev=0.0016106433  ; 88.1 km at 78N
;targetLev=0.001121129  ; 89.5 km at 78N
targetLev=0.005842710  ; 83 km at 78N
;targetLev=0.00395005   ; 84.3 km at 82N

caseid="125b36Pb02"
;caseid="125b36eg"

coreDistEnsembleMeso(targetLat, targetLev, caseid)
