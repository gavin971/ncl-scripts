load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;path="$MODEL_PATH/ccsm/case/"
path="$MODEL_PATH/"

cntlCase="b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678"
exptCase="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"

;yearstring="2014-2018"
;yearstring="2017"
;yearstring="2018"
yearstring="2014"

doDJF=False

if (doDJF) then
  cntlFilepathDJF=path+cntlCase+"/proc/"+cntlCase+".cam2.h0."+yearstring+"DJF.nc"
  exptFilepathDJF=path+exptCase+"/proc/"+exptCase+".cam2.h0."+yearstring+"DJF.nc"

  cntlFileDJF=addfile(cntlFilepathDJF,"r")
  exptFileDJF=addfile(exptFilepathDJF,"r")

  cPRECTDJF=cntlFileDJF->PRECT
  ePRECTDJF=exptFileDJF->PRECT
  dPRECTDJF=ePRECTDJF ; get metadata
  dPRECTDJF=dPRECTDJF-cPRECTDJF
  dPRECTDJFpct=dPRECTDJF ; get metadata
  dPRECTDJFpct=100.*dPRECTDJFpct/cPRECTDJF
  ; Convert m/s to mm/day
  dPRECTDJF=dPRECTDJF*1000.*24.*3600.
end if

cntlFilepathJJA=path+cntlCase+"/proc/"+cntlCase+".cam2.h0."+yearstring+"JJA.nc"
exptFilepathJJA=path+exptCase+"/proc/"+exptCase+".cam2.h0."+yearstring+"JJA.nc"

cntlFileJJA=addfile(cntlFilepathJJA,"r")
exptFileJJA=addfile(exptFilepathJJA,"r")

cPRECTJJA=cntlFileJJA->PRECT
ePRECTJJA=exptFileJJA->PRECT
dPRECTJJA=ePRECTJJA ; get metadata
dPRECTJJA=dPRECTJJA-cPRECTJJA
dPRECTJJApct=dPRECTJJA ; get metadata
dPRECTJJApct=100.*dPRECTJJApct/cPRECTJJA

; Convert m/s to mm/day
dPRECTJJA=dPRECTJJA*1000.*24.*3600.

;************************************************
; create plot
;************************************************
;wks = gsn_open_wks("eps", "$IMG_PATH/nuclearTSdiff")                ; open an eps file
wks = gsn_open_wks("pdf", "$IMG_PATH/nuclearPRECTdiff"+yearstring)                ; open a pdf file

nPanel=2
plot = new(2, graphic)

res                       = True     ; plot mods desired
res@gsnDraw   = False
res@gsnFrame  = False
;res@vpWidthF  = 1.0		; set width of plot 
;res@vpHeightF = 1.0		; set height of plot

res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = True    ; turn of contour lines
res@cnLevelSelectionMode  = "ExplicitLevels"
res@cnLevels = (/-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2/)
;res@cnLevelSelectionMode  = "ManualLevels"
;res@cnLevelSpacingF = 0.5
;res@cnMinLevelValF  = -3.5
;res@cnMaxLevelValF  = 3.0

res@gsnSpreadColors       = False     ; use full range of color map

res@lbLabelStride         = 1
res@lbOrientation         = "Vertical"

res@pmTickMarkDisplayMode = "Always"; use NCL default lat/lon labels

;  res@gsnAddCyclic          = False    ; data already has cyclic point
				     ; this must also be set for any zoom

; note that the gsn_csm_*map_ce templates automatically set 
; res@mpLimitMode="LatLon" for you. If you are plotting a different projection,
; you may have to set this resource.


;gsn_define_colormap(wks,"BlWhRe")                  ; change colormap 
; Change color map using named colors
  colors = (/"white","black",\
             "sienna4","darkgoldenrod4","darkgoldenrod3","darkgoldenrod1", \
             "white", \
             "chartreuse1","chartreuse2","chartreuse3","chartreuse4"/)

  gsn_define_colormap(wks,colors)                  ; change colormap 


res@gsnLeftString        = " "
res@gsnRightString       = " "
res@gsnCenterString      = " "
res@tmXBLabelsOn         = True
res@lbLabelBarOn         = True
res@lbOrientation        = "Vertical"
;res@lbTopMarginF         = 0.0
;res@lbBottomMarginF      = -2.8
;res@lbRightMarginF         = 0.3
res@tmXBMajorLengthF     = 0.0
plot(0) = gsn_csm_contour_map_ce(wks,dPRECTJJA(0,:,:), res)

;delete(res@cnLevels)
res@cnLevels = (/-80, -60,-40,-20,20,40,60,80/)
res@lbLabelBarOn         = True
;res@lbTopMarginF         = 0.0
;res@lbBottomMarginF      = -2.8
res@tmXBLabelsOn         = False
plot(1) = gsn_csm_contour_map_ce(wks,dPRECTJJApct(0,:,:), res)

pres                        = True
pres@gsnMaximize            = True
pres@gsnPaperOrientation    = "landscape"
gsn_panel(wks, plot, (/ 2, 1 /), pres)


;res@cnLevels = (/-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2/)
;res@tmXBLabelsOn         = False
;res@lbLabelBarOn         = True
;res@lbOrientation        = "Vertical"
;res@lbLabelFontHeightF   = 0.01
;plot(0) = gsn_csm_contour_map_ce(wks,dPRECTDJF(0,:,:), res)
;res@cnLevels = (/-80, -60,-40,-20,20,40,60,80/)
;res@tmXBLabelsOn         = False
;res@lbLabelBarOn         = True
;res@lbOrientation        = "Vertical"
;res@lbLabelFontHeightF   = 0.01
;plot(1) = gsn_csm_contour_map_ce(wks,dPRECTDJFpct(0,:,:), res)
;gsn_panel(wks, plot, (/ 2, 1 /), pres)

