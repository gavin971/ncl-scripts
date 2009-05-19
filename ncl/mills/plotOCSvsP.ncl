;*************************************************
; plotOCSvsN2O.ncl
; by Mike Mills
; May 15, 2009
;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "tests/common.ncl"
load "tests/getHistoryPath.ncl"
load "tests/getDataPath.ncl"
;*************************************************

begin
  dataPath=getDataPath()
  ATLAS08_20N = readAsciiTable(dataPath+"ATMOS/SPARC_Fig2.9/8-20N.txt", 4,"float",1)
  ATLAS20_30N = readAsciiTable(dataPath+"ATMOS/SPARC_Fig2.9/20-30N.txt",2,"float",1)
  ATLAS30_40N = readAsciiTable(dataPath+"ATMOS/SPARC_Fig2.9/30-40N.txt",2,"float",1)
  ATLAS40_49N = readAsciiTable(dataPath+"ATMOS/SPARC_Fig2.9/40-49N.txt",2,"float",1)
  
  historyRun="125smin"
  historyPattern=getHistoryPath(historyRun)+"daily/*.cam2.h3.1997-11*.nc"
  historyNames = systemfunc("ls " + historyPattern)
  print(historyNames)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  OCS=addfiles_GetVar(history, historyNames, "OCS")
; OCS[time | 55] x [lev | 66] x [lat | 46] x [lon | 72]
  time = addfiles_GetVar(history, historyNames, "time")
;  print(time)
  lev = history[0]->lev
  lat = history[0]->lat
  lon = history[0]->lon
  nlev=dimsizes(lev)

  OCSmean_20_30N=new(nlev,float)
  OCSmean_20_30N!0="lev"
  OCSmean_20_30N&lev=lev
  
  OCSmean_8_20N=OCSmean_20_30N
  OCSmean_30_40N=OCSmean_20_30N
  OCSmean_40_49N=OCSmean_20_30N
  
  OCSerr_8_20N=OCSmean_20_30N
  OCSerr_20_30N=OCSmean_20_30N
  OCSerr_30_40N=OCSmean_20_30N
  OCSerr_40_49N=OCSmean_20_30N
  
  ; get mean and standard deviation at each level
  ; for November, in each lat band, all longitudes
  do i= 0, nlev-1
    OCSmean_8_20N(i)=avg(OCS(:,i,{8:20},:))
    OCSerr_8_20N(i)=stddev(OCS(:,i,{8:20},:))
    OCSmean_20_30N(i)=avg(OCS(:,i,{20:30},:))
    OCSerr_20_30N(i)=stddev(OCS(:,i,{20:30},:))
    OCSmean_30_40N(i)=avg(OCS(:,i,{30:40},:))
    OCSerr_30_40N(i)=stddev(OCS(:,i,{30:40},:))
    OCSmean_40_49N(i)=avg(OCS(:,i,{40:49},:))
    OCSerr_40_49N(i)=stddev(OCS(:,i,{40:49},:))
  end do
  OCSmean_8_20N=OCSmean_8_20N*1e12 ; pptv
  OCSerr_8_20N=OCSerr_8_20N*1e12 ; pptv
  OCSmean_20_30N=OCSmean_20_30N*1e12 ; pptv
  OCSerr_20_30N=OCSerr_20_30N*1e12 ; pptv
  OCSmean_30_40N=OCSmean_30_40N*1e12 ; pptv
  OCSerr_30_40N=OCSerr_30_40N*1e12 ; pptv
  OCSmean_40_49N=OCSmean_40_49N*1e12 ; pptv
  OCSerr_40_49N=OCSerr_40_49N*1e12 ; pptv

;************************************************
; plotting parameters
;************************************************
  wks   = gsn_open_wks ("pdf","OCSvsP_"+historyRun)                 ; open workstation

  res                   = True                     ; plot mods desired
  res@gsnDraw      = False                        ; don't draw yet
  res@gsnFrame     = False                        ; don't advance frame yet
;  res@tiMainString      = "Scatter Plot"           ; add title
  res@xyMarkLineMode   = "Markers"                ; choose which have markers
  res@xyMarker         =  1                      ; choose type of marker (dots)
;  res@xyMarkerSizeF     = 0.01                    ; Marker size (default 0.01)
  res@trXMinF           = 0.0
  res@trXMaxF           = 550.0
  res@trYMinF           = 5.0
  res@trYMaxF           = 500.0
  res@trYLog            = True                       ; log Y-axis
  res@trYReverse        = True                       ; reverse Y-axis
  res@tiXAxisString     = "OCS (pptv)"
  res@tiYAxisString     = "Pressure (hPa)"

  res@tmLabelAutoStride = True                    ; nice tick mark labels

; Make custom square markers:
;  res@xyMarkers     = NhlNewMarker(wks, "^", 19, 0.0, 0.0, 1.5, 1.5, 0.0)
;  res@xyMarkers     = 16
  res@xyMarkLineMode   = "Lines"                ; choose lines
  res@xyMonoLineColor     = False
  res@xyLineColors     = (/ "Blue", "Green", "Red", "Black", "Blue", "Blue", "Green", "Green", "Red", "Red", "Black", "Black"/)                  ; Line color
  res@xyLineThicknesses = (/ 3.0, 3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 /)
  res@xyMonoLineThickness = False
  
  data2      = new((/12,nlev/),float)

  data2(0,:) = OCSmean_8_20N
  data2(1,:) = OCSmean_20_30N
  data2(2,:) = OCSmean_30_40N
  data2(3,:) = OCSmean_40_49N
  data2(4,:) = OCSmean_8_20N-OCSerr_8_20N
  data2(5,:) = OCSmean_8_20N+OCSerr_8_20N
  data2(6,:) = OCSmean_20_30N-OCSerr_20_30N
  data2(7,:) = OCSmean_20_30N+OCSerr_20_30N
  data2(8,:) = OCSmean_30_40N-OCSerr_30_40N
  data2(9,:) = OCSmean_30_40N+OCSerr_30_40N
  data2(10,:) = OCSmean_40_49N-OCSerr_40_49N
  data2(11,:) = OCSmean_40_49N+OCSerr_40_49N

  data2!1="lev"
  data2&lev=lev
  
;************************************************
; Add altitude tick labels on the right Y-axis:  
;************************************************
  add_hgt=False
  set_pres_hgt_axes(lev,res,add_hgt)
  
  res@xyDashPattern = 0               ; Make curves all solid
 
  plot  = gsn_csm_xy (wks,data2,lev,res) ; create plot

;************************************************
; Add a right Y axis "height" label.
;************************************************
  if(add_hgt)
    rightaxis_string = create "right_axis" textItemClass wks
      "txString"      : "Height (km)"
      "txFontHeightF" : 0.025
      "txAngleF"      : 90.
    end create

    anno = NhlAddAnnotation(plot,rightaxis_string)

    setvalues anno
      "amZone"          : 3      ; Just outside plot area
      "amJust"          : "centercenter"
      "amSide"          : "right" 
      "amParallelPosF"  : 0.5
      "amOrthogonalPosF": 0.03
      "amResizeNotify"  : True     ; Resize if plot resized.
    end setvalues
  end if

;************************************************
; add observations
;************************************************  

  res@xyMonoLineColor     = True
  res@xyDashPattern = 1               ; Make curves dashed
  
  res@xyLineColor     = "Blue"                  ; Line color
  overlay(plot,gsn_csm_xy (wks,ATLAS08_20N(:,0),ATLAS08_20N(:,1),res))
;  print(ATLAS08_20N(:,0)+" | "+ATLAS08_20N(:,1))

  res@xyLineColor     = "Green"                 ; Line color
  overlay(plot,gsn_csm_xy (wks,ATLAS20_30N(:,0),ATLAS20_30N(:,1),res))
;  print(ATLAS20_30N(:,0)+" | "+ATLAS20_30N(:,1))

  res@xyLineColor     = "Red"                   ; Line color
  overlay(plot,gsn_csm_xy (wks,ATLAS30_40N(:,0),ATLAS30_40N(:,1),res))
;  print(ATLAS30_40N(:,0)+" | "+ATLAS30_40N(:,1))

  res@xyLineColor     = "Black"                 ; Line color
  overlay(plot,gsn_csm_xy (wks,ATLAS40_49N(:,0),ATLAS40_49N(:,1),res))
;  print(ATLAS40_49N(:,0)+" | "+ATLAS40_49N(:,1))
  
  legend = create "Legend" legendClass wks 
    "vpXF"                     : 0.55                   ; coordinate of left side of legend
    "vpYF"                     : 0.8                    ; coordinate of top of legend
    "vpWidthF"                 : 0.15                     ; width
    "vpHeightF"                : 0.2                     ; height
    "lgPerimOn"                : False                   ; no perimeter
    "lgItemCount"              : 4                      ; how many
    "lgLabelStrings"           : (/ "8-20~S~o~N~N", "20-30~S~o~N~N", "30-40~S~o~N~N", "40-49~S~o~N~S"/)
    "lgLabelsOn"               : True
    "lgMonoItemType"           : True
    "lgItemType"               : "Lines"
    "lgMonoLineColor"          : False
    "lgLineColors"             : (/ "blue", "green", "red", "black" /)   ; markers colors
    "lgMonoDashIndex"          : True ; All solid lines
    "lgDashIndex"              : "SolidLine"
    "lgLineThicknessF"         : 3.0
    "lgMonoLineLabelFontColor" : True                    ; one label color
    "lgAutoManage"             : False
    "lgItemPlacement"          : "ExplicitPlacement"
    "lgItemPositions"          : (/ 8.0, 8.5, 9.0, 9.5 /)
    "lgItemOrder"              : (/ 3, 2, 1, 0 /)
    "lgLabelPosition"          : "Right"
    "lgLabelJust"              : "CenterLeft"
  end create
    
; these two step are required!
  draw(plot)                     ; note we are drawing the first one!
  draw(legend)                   ; add legend
  frame(wks)

end
