load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

path="$MODEL_PATH/ccsm/case/"

cntlCase="b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678"
exptCase="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"

cntlFilepathDJF=path+cntlCase+"/proc/"+cntlCase+".cam2.h0.2014-2018DJF.nc"
exptFilepathDJF=path+exptCase+"/proc/"+exptCase+".cam2.h0.2014-2018DJF.nc"

cntlFileDJF=addfile(cntlFilepathDJF,"r")
exptFileDJF=addfile(exptFilepathDJF,"r")

cTSDJF=cntlFileDJF->TS
eTSDJF=exptFileDJF->TS
dTSDJF=eTSDJF ; get metadata
dTSDJF=dTSDJF-cTSDJF

cntlFilepathJJA=path+cntlCase+"/proc/"+cntlCase+".cam2.h0.2014-2018JJA.nc"
exptFilepathJJA=path+exptCase+"/proc/"+exptCase+".cam2.h0.2014-2018JJA.nc"

cntlFileJJA=addfile(cntlFilepathJJA,"r")
exptFileJJA=addfile(exptFilepathJJA,"r")

cTSJJA=cntlFileJJA->TS
eTSJJA=exptFileJJA->TS
dTSJJA=eTSJJA ; get metadata
dTSJJA=dTSJJA-cTSJJA


;************************************************
; create plot
;************************************************
;wks = gsn_open_wks("eps", "$IMG_PATH/nuclearTSdiff")                ; open an eps file
wks = gsn_open_wks("pdf", "$IMG_PATH/nuclearTSdiff")                ; open a pdf file

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
res@cnLevels = (/-8,-6,-5,-4, -3, -2.5,-2, -1.5, -1, -0.5, 0.5, 1/)

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
  colors = (/"white","black","darkorchid","orchid","blue","dodgerblue3", \
             "dodgerblue2","dodgerblue1","cyan4","cyan3","cyan2","cyan1","white", \
             "yellow","orange"/)

  gsn_define_colormap(wks,colors)                  ; change colormap 


res@gsnLeftString        = " "
res@gsnRightString       = " "
res@gsnCenterString      = " "
res@tmXBLabelsOn         = True
res@lbLabelBarOn         = True
res@lbTopMarginF         = 0.0
res@lbBottomMarginF      = -2.8
;res@lbRightMarginF         = 0.3
res@tmXBMajorLengthF     = 0.0
plot(0) = gsn_csm_contour_map_ce(wks,dTSJJA(0,:,:), res)

res@tmXBLabelsOn         = False
res@lbLabelBarOn         = True
res@lbOrientation        = "Vertical"
res@lbTopMarginF         = 5
res@lbBottomMarginF      = -7
res@lbLabelFontHeightF   = 0.01
;delete(res@tmXBMajorLengthF)
plot(1) = gsn_csm_contour_map_ce(wks,dTSDJF(0,:,:), res)

pres                        = True
pres@gsnMaximize            = True
pres@gsnPaperOrientation    = "landscape"
gsn_panel(wks, plot, (/ 2, 1 /), pres)


