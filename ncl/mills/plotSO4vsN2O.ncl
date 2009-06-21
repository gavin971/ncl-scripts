;*************************************************
; plotSO4vsN2O.ncl
; by Mike Mills
; June 10, 2009
;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "tests/getHistoryPath.ncl"
load "tests/getDataPath.ncl"
;*************************************************

begin
  dataPath=getDataPath()
  wilson = readAsciiTable(dataPath+"Chuck Wilson/WilsonSO4vsN2O.txt",4,"float",1)

  pattern="h0"
  subDir="allmonthly"

  historyRun="nogeo400SO2"
  historyPath = getHistoryPath(historyRun)
  historyPattern=historyPath+subDir+"/*.cam2."+pattern+"*.nc"
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  print("historyNames:"+dimsizes(historyNames))
  N2O = addfiles_GetVar(history, historyNames, "N2O")
  PURSULVMR = addfiles_GetVar(history, historyNames, "PURSULVMR")
  MIXSULVMR = addfiles_GetVar(history, historyNames, "MIXSULVMR")
  SULVMR=PURSULVMR ; get metadata
  SULVMR=SULVMR+MIXSULVMR
  
;  historyRun="nogeoPurSpin1"
;  historyPath = getHistoryPath(historyRun)
;  historyPattern=historyPath+subDir+"/*.cam2."+pattern+"*.nc"
;  historyNames2 = systemfunc("ls " + historyPattern)
;  history2 = addfiles(historyNames2, "r")
;  ListSetType(history2, "cat")
;  print("historyNames2:"+dimsizes(historyNames2))
;  N2O_2 = addfiles_GetVar(history2, historyNames2, "N2O")
;  PURSULVMR_2 = addfiles_GetVar(history2, historyNames2, "PURSULVMR")
;  MIXSULVMR_2 = addfiles_GetVar(history2, historyNames2, "MIXSULVMR")
;  SULVMR_2=PURSULVMR_2 ; get metadata
;  SULVMR_2=SULVMR_2+MIXSULVMR_2
  
; [time | 13] x [lev | 125] x [lat | 46] x [lon | 72]

;************************************************
; plotting parameters
;************************************************
  wks   = gsn_open_wks ("pdf","SO4vsN2O")                 ; open workstation

  res                   = True                     ; plot mods desired
  res@gsnDraw      = False                        ; don't draw yet
  res@gsnFrame     = False                        ; don't advance frame yet
;  res@tiMainString      = "Scatter Plot"           ; add title
  res@xyMarkLineMode   = "Markers"                ; choose which have markers
  res@xyMarker         =  1                      ; choose type of marker (dots)
;  res@xyMarkerSizeF     = 0.01                    ; Marker size (default 0.01)
  res@trXMinF           = 40.0
  res@trXMaxF           = 340.0
  res@trYMinF           = 0.02
  res@trYMaxF           = 2.0
  res@tiXAxisString     = "N~B~2~N~O (ppbv)"
  res@tiYAxisString     = "SO~B~4~N~ (ppbv)"
  res@trYLog            = True

  res@tmLabelAutoStride = True                    ; nice tick mark labels

; Make custom square markers:
;  res@xyMarkers     = NhlNewMarker(wks, "^", 19, 0.0, 0.0, 1.5, 1.5, 0.0)
;  res@xyMarkers     = 16
  res@xyMarkLineMode   = "Lines"                ; choose lines
  res@xyMonoLineColor     = True
  res@xyLineColor     = "Black"                    ; Line color
  
  plot  = gsn_csm_xy (wks,N2O(0,{3:1000},0,:)*1e9,SULVMR(0,{3:1000},0,:)*1e9,res) ; create plot
  i = 8
    do j = 0,45
      print(i+" "+j)
      overlay(plot,gsn_csm_xy (wks,N2O(i,{3:1000},j,:)*1e9,SULVMR(i,{3:1000},j,:)*1e9,res))
    end do
;  end do
  
;  res@xyLineColor     = "Chartreuse4"                    ; Line color
;  i = 15
;  do j = 0,45
;    print(i+" "+j)  
;    overlay(plot,gsn_csm_xy (wks,N2O_2(i,{3:1000},j,:)*1e9,SULVMR_2(i,{3:1000},j,:)*1e9,res))
;  end do
  
;==========================
; Wilson data
;==========================  
  res@xyMarkLineMode   = "MarkLines"             ; choose markers
  res@xyMarkerColor    = "red"                  ; Marker color
  res@xyLineColor      = "red"                    ; Line color
  res@xyMarkers        =  16                    ; choose type of marker (closed circles)
  res@xyLineThicknessF = 4
  overlay(plot,gsn_csm_xy (wks,wilson(:,0),wilson(:,1),res))

;==========================
; Error bars (Wilson data)
;==========================  
  ndims=dimsizes(wilson)
; gsn_add* templates are functions that we set to dummy values. Since
; we are going to draw numerous error bars, we create two arrays to
; hold the dummy values.
  error_bar = new(ndims(0),graphic)
;  print(ndims)
  
  polyres                   = True                       ; marker resources
  polyres@gsLineThicknessF  = 4
  polyres@gsLineColor      = "red"                    ; Line color

  do t=0,ndims(0)-1
;    print(t+" "+wilson(t,:))
    xarr=(/wilson(t,0), wilson(t,0)/)
    yarr=(/wilson(t,1) + wilson(t,2), wilson(t,1) + wilson(t,3)/)
;    print(error_bar(t))
    error_bar(t) = gsn_add_polyline (wks, plot, xarr, yarr, polyres)
  end do
  
; these two step are required!
  draw(plot)                     ; note we are drawing the first one!
;  draw(legend)                   ; add legend
  frame(wks)

end