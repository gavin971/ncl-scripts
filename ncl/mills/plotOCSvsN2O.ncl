;*************************************************
; plotOCSvsN2O.ncl
; by Mike Mills
; May 15, 2009
;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "tests/getHistoryPath.ncl"
;*************************************************

begin
  Manitoba19960724 = readAsciiTable("~/Data/Geoff Toon OCSvN2O/96.07.24Manitoba.txt",5,"float",2)
  NewMexico19960928 = readAsciiTable("~/Data/Geoff Toon OCSvN2O/96.09.28NewMexico.txt",5,"float",2)
  Fairbanks19970508 = readAsciiTable("~/Data/Geoff Toon OCSvN2O/97.05.08Fairbanks.txt",5,"float",2)
  Fairbanks19970708a = readAsciiTable("~/Data/Geoff Toon OCSvN2O/97.07.08FairbanksAscent.txt",5,"float",2)
  Fairbanks19970708d = readAsciiTable("~/Data/Geoff Toon OCSvN2O/97.07.08FairbanksDescent.txt",5,"float",2)

  filename=getHistoryPath("125smin")+"monthly/all.125smin.cam2.h0.1998-04.nc"
  file1=addfile(filename,"r")
  N2O=file1->N2O
  OCS=file1->OCS

;************************************************
; plotting parameters
;************************************************
  wks   = gsn_open_wks ("ps","xy")                 ; open workstation

  res                   = True                     ; plot mods desired
  res@gsnDraw      = False                        ; don't draw yet
  res@gsnFrame     = False                        ; don't advance frame yet
;  res@tiMainString      = "Scatter Plot"           ; add title
  res@xyMarkLineModes   = "Markers"                ; choose which have markers
  res@xyMarkers         =  1                      ; choose type of marker (dots)
  res@xyMarkerColor     = "red"                    ; Marker color
  res@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
  res@trXMinF           = 0.0
  res@trXMaxF           = 350.0
  res@trYMinF           = 0.0
  res@trYMaxF           = 550.0
  res@tiXAxisString     = "N~B~2~N~O (ppbv)"
  res@tiYAxisString     = "OCS (pptv)"

  res@tmLabelAutoStride = True                    ; nice tick mark labels

  i = 30
  plot  = gsn_csm_xy (wks,N2O(0,:,i,:)*1e9,OCS(0,:,i,:)*1e12,res) ; create plot
  do i = 31,40
    overlay(plot,gsn_csm_xy (wks,N2O(0,:,i,:)*1e9,OCS(0,:,i,:)*1e12,res))
  end do
  
  do i = 5,15
    overlay(plot,gsn_csm_xy (wks,N2O(0,:,i,:)*1e9,OCS(0,:,i,:)*1e12,res))
  end do
  
  res@xyMarkers         =  16                      ; choose type of marker (closed circles)
  overlay(plot,gsn_csm_xy (wks,Manitoba19960724(:,3)*1e9,Manitoba19960724(:,1)*1e12,res))

  res@xyMarkerColor     = "blue"                  ; Marker color
  overlay(plot,gsn_csm_xy (wks,NewMexico19960928(:,3)*1e9,NewMexico19960928(:,1)*1e12,res))

  res@xyMarkerColor     = "green"                 ; Marker color
  overlay(plot,gsn_csm_xy (wks,Fairbanks19970508(:,3)*1e9,Fairbanks19970508(:,1)*1e12,res))

  res@xyMarkerColor     = "black"                 ; Marker color
  overlay(plot,gsn_csm_xy (wks,Fairbanks19970708a(:,3)*1e9,Fairbanks19970708a(:,1)*1e12,res))

  res@xyMarkerColor     = "purple"                ; Marker color
  overlay(plot,gsn_csm_xy (wks,Fairbanks19970708d(:,3)*1e9,Fairbanks19970708d(:,1)*1e12,res))
    
; these two step are required!
  draw(plot)                     ; note we are drawing the first one!
  frame(wks)

end
