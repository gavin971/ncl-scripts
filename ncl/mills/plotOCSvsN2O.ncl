;*************************************************
; plotOCSvsN2O.ncl
; by Mike Mills
; May 15, 2009
;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "tests/getHistoryPath.ncl"
load "tests/getDataPath.ncl"
;*************************************************

begin
  dataPath=getDataPath()
  Manitoba19960724 = readAsciiTable(dataPath+"Geoff Toon OCSvN2O/96.07.24Manitoba.txt",5,"float",2)
  NewMexico19960928 = readAsciiTable(dataPath+"Geoff Toon OCSvN2O/96.09.28NewMexico.txt",5,"float",2)
  Fairbanks19970508 = readAsciiTable(dataPath+"Geoff Toon OCSvN2O/97.05.08Fairbanks.txt",5,"float",2)
  Fairbanks19970708a = readAsciiTable(dataPath+"Geoff Toon OCSvN2O/97.07.08FairbanksAscent.txt",5,"float",2)
  Fairbanks19970708d = readAsciiTable(dataPath+"Geoff Toon OCSvN2O/97.07.08FairbanksDescent.txt",5,"float",2)

  filename=getHistoryPath("125smin")+"monthly/all.125smin.cam2.h0.1997-07.nc"
  file1=addfile(filename,"r")
  N2O=file1->N2O ; [time | 1] x [lev | 125] x [lat | 46] x [lon | 72]
  OCS=file1->OCS

;************************************************
; plotting parameters
;************************************************
  wks   = gsn_open_wks ("pdf","OCSvsN2O")                 ; open workstation

  res                   = True                     ; plot mods desired
  res@gsnDraw      = False                        ; don't draw yet
  res@gsnFrame     = False                        ; don't advance frame yet
;  res@tiMainString      = "Scatter Plot"           ; add title
  res@xyMarkLineMode   = "Markers"                ; choose which have markers
  res@xyMarker         =  1                      ; choose type of marker (dots)
;  res@xyMarkerSizeF     = 0.01                    ; Marker size (default 0.01)
  res@trXMinF           = 0.0
  res@trXMaxF           = 350.0
  res@trYMinF           = 0.0
  res@trYMaxF           = 550.0
  res@tiXAxisString     = "N~B~2~N~O (ppbv)"
  res@tiYAxisString     = "OCS (pptv)"

  res@tmLabelAutoStride = True                    ; nice tick mark labels

; Make custom square markers:
;  res@xyMarkers     = NhlNewMarker(wks, "^", 19, 0.0, 0.0, 1.5, 1.5, 0.0)
;  res@xyMarkers     = 16
  res@xyMarkLineMode   = "Lines"                ; choose lines
  res@xyMonoLineColor     = True
  res@xyLineColor     = "Black"                    ; Line color
  
  i = 30
  plot  = gsn_csm_xy (wks,N2O(0,:,i,:)*1e9,OCS(0,:,i,:)*1e12,res) ; create plot
  do i = 31,40
    overlay(plot,gsn_csm_xy (wks,N2O(0,:,i,:)*1e9,OCS(0,:,i,:)*1e12,res))
  end do
  
  res@xyLineColor     = "Chartreuse4"                    ; Line color  
  do i = 5,15
    overlay(plot,gsn_csm_xy (wks,N2O(0,:,i,:)*1e9,OCS(0,:,i,:)*1e12,res))
  end do
  
  res@xyMarkLineMode   = "Markers"                ; choose markers
  res@xyMarkerColor     = "red"                    ; Marker color
  res@xyMarkers         =  16                      ; choose type of marker (closed circles)
  overlay(plot,gsn_csm_xy (wks,Manitoba19960724(:,3)*1e9,Manitoba19960724(:,1)*1e12,res))

  res@xyMarkerColor     = "blue"                  ; Marker color
  overlay(plot,gsn_csm_xy (wks,NewMexico19960928(:,3)*1e9,NewMexico19960928(:,1)*1e12,res))

  res@xyMarkerColor     = "green"                 ; Marker color
  overlay(plot,gsn_csm_xy (wks,Fairbanks19970508(:,3)*1e9,Fairbanks19970508(:,1)*1e12,res))

  res@xyMarkerColor     = "purple"                 ; Marker color
  overlay(plot,gsn_csm_xy (wks,Fairbanks19970708a(:,3)*1e9,Fairbanks19970708a(:,1)*1e12,res))

  overlay(plot,gsn_csm_xy (wks,Fairbanks19970708d(:,3)*1e9,Fairbanks19970708d(:,1)*1e12,res))

  legend = create "Legend" legendClass wks 
    "vpXF"                     : 0.21                   ; coordinate of left side of legend
    "vpYF"                     : 0.8                    ; coordinate of top of legend
    "vpWidthF"                 : 0.15                     ; width
    "vpHeightF"                : 0.2                     ; height
    "lgPerimOn"                : False                   ; no perimeter
    "lgItemCount"              : 6                      ; how many
    "lgLabelStrings"           : (/ "Manitoba 7/24/1996", "New Mexico 9/28/1996", "Fairbanks 5/8/1997", "Fairbanks 7/8/1997", "WACCM, July 30-70~S~o~N~N", "WACCM, July 30-70~S~o~N~S"/)
    "lgLabelsOn"               : True
    "lgMonoItemType"           : False
    "lgItemTypes"               : (/ "Markers", "Markers", "Markers", "Markers", "Lines", "Lines" /)
    "lgMonoMarkerIndex"        : True
    "lgMarkerIndex"            : 16                     ; markers
    "lgMonoMarkerColor"        : False
    "lgMarkerColors"           : (/ "red", "blue", "green", "purple", "black", "Chartreuse4" /)   ; markers colors
    "lgMonoLineColor"          : False
    "lgLineColors"             : (/ "red", "blue", "green", "purple", "black", "Chartreuse4" /)   ; markers colors
    "lgMonoDashIndex"          : True ; All solid lines
    "lgDashIndex"              : "SolidLine"
    "lgLineThicknessF"         : 3.0
    "lgMonoLineLabelFontColor" : True                    ; one label color
    "lgAutoManage"             : False
    "lgItemPlacement"          : "ExplicitPlacement"
    "lgItemPositions"          : (/ 7.0, 7.5, 8.0, 8.5, 9.0, 9.5 /)
    "lgItemOrder"              : (/ 5, 4, 3, 2, 1, 0 /)
    "lgLabelPosition"          : "Right"
    "lgLabelJust"              : "CenterLeft"
  end create
    
; these two step are required!
  draw(plot)                     ; note we are drawing the first one!
  draw(legend)                   ; add legend
  frame(wks)

end
