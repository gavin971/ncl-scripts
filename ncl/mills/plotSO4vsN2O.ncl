;*********************************************************
; plotSO4vsN2O.ncl
; by Mike Mills
; June 10, 2009
;*********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "ncl-scripts/getHistoryPath.ncl"
load "ncl-scripts/getDataPath.ncl"

;*********************************************************
; Run allSaveSulfVMR.ncl prior to this script to generate
; the analysis files read in by this script. 
;*********************************************************

undef("plotSO4vsN2O")
procedure plotSO4vsN2O (historyRun1:string, historyRun2:string)
begin
  dataPath = getDataPath()
  dataFile=dataPath+"Chuck Wilson/WilsonSO4vsN2O.txt"
  print(dataFile)
  wilson = readAsciiTable(dataFile,4,"float",1)

  subDir="allmonthly"

  historyPath = "analysis/"+historyRun1+"/"
  historyFile = addfile(historyPath+"sulfvmr."+historyRun1+".nc","r")
  N2O = historyFile->N2O
  PURSULVMR = historyFile->PURSULVMR
  MIXSULVMR = historyFile->MIXSULVMR
  MIXSULSHELLVMR = historyFile->MIXSULSHELLVMR
  SULVMR=PURSULVMR ; get metadata
  SULVMR=SULVMR+MIXSULSHELLVMR
;  SULVMR=SULVMR+MIXSULVMR
  SULVMR@_FillValue=0.0
  nTimes=dimsizes(historyFile->time)
  
  historyPath = "analysis/"+historyRun2+"/"
  historyFile2 = addfile(historyPath+"sulfvmr."+historyRun2+".nc","r")
  N2O_2 = historyFile2->N2O
  PURSULVMR_2 = historyFile2->PURSULVMR
  MIXSULVMR_2 = historyFile2->MIXSULVMR
  MIXSULSHELLVMR_2 = historyFile2->MIXSULSHELLVMR
  SULVMR_2=PURSULVMR_2 ; get metadata
  SULVMR_2=SULVMR_2+MIXSULSHELLVMR_2
;  SULVMR_2=SULVMR_2+MIXSULVMR_2
  SULVMR_2@_FillValue=0.0
  nTimes2=dimsizes(historyFile2->time)
  
; [time | 13] x [lev | 125] x [lat | 46] x [lon | 72]

;************************************************
; plotting parameters
;************************************************
  wks   = gsn_open_wks ("pdf","images/SO4vsN2O."+historyRun1+"."+historyRun2)                 ; open workstation

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
  res@trYMaxF           = 20.0
  res@tiXAxisString     = "N~B~2~N~O (ppbv)"
  res@tiYAxisString     = "SO~B~4~N~ (ppbv)"
  res@trYLog            = True

  res@tmLabelAutoStride = True                    ; nice tick mark labels

; Make custom square markers:
;  res@xyMarkers     = NhlNewMarker(wks, "^", 19, 0.0, 0.0, 1.5, 1.5, 0.0)
;  res@xyMarkers     = 16
  res@xyMarkLineMode   = "Lines"                ; choose lines
  res@xyMonoLineColor     = True
  res@xyLineColor     = "red"                    ; Line color
  
  plot  = gsn_csm_xy (wks,N2O(0,{40:250},0,:)*1e9,SULVMR(0,{40:250},0,:)*1e9,res) ; create plot
  do i = 0, nTimes-1
    print((i+1)+" of "+nTimes)
    do j = 21,45  ; 5S - 90N
;      print(i+" "+j)
      overlay(plot,gsn_csm_xy (wks,N2O(i,{40:250},j,:)*1e9,SULVMR(i,{40:250},j,:)*1e9,res))
    end do
  end do
  
  res@xyLineColor     = "black"                    ; Line color
  do i = 0, nTimes2-1
    print((i+1)+" of "+nTimes2)
    do j = 21,45 ; 5S - 90N
;      print(i+" "+j)  
      overlay(plot,gsn_csm_xy (wks,N2O_2(i,{40:250},j,:)*1e9,SULVMR_2(i,{40:250},j,:)*1e9,res))
    end do
  end do
  
  print("plotting Wilson data...")
  
;==========================
; Wilson data
;==========================  
  res@xyMarkLineMode   = "MarkLines"             ; choose markers
  res@xyMarkerColor    = "blue"                  ; Marker color
  res@xyLineColor      = "blue"                    ; Line color
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
  polyres@gsLineColor      = "blue"                    ; Line color

  do t=0,ndims(0)-1
;    print(t+" "+wilson(t,:))
    xarr=(/wilson(t,0), wilson(t,0)/)
    yarr=(/wilson(t,1) + wilson(t,2), wilson(t,1) + wilson(t,3)/)
;    print(error_bar(t))
    error_bar(t) = gsn_add_polyline (wks, plot, xarr, yarr, polyres)
  end do
  
; these two step are required!
  print("drawing plot...")
  draw(plot)                     ; note we are drawing the first one!
;  draw(legend)                   ; add legend
  print("frame...")
  frame(wks)

end

;*************************************************

plotSO4vsN2O ("125b36eg", "125b36Pb02")
