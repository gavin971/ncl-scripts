load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

path="$MODEL_PATH/ccsm/case/"

cntlCase="b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678"
exptCase="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"


;cntlFilepath=path+cntlCase+"/proc/"+cntlCase+".cam.growseas.2014-2018avg.nc"
;exptFilepath=path+exptCase+"/proc/"+exptCase+".cam.growseas.2014-2018avg.nc"

cntlFilepath=path+cntlCase+"/proc/"+cntlCase+".cam.growseas.2014-2021avg.nc"
exptFilepath=path+exptCase+"/proc/"+exptCase+".cam.growseas.2014-2021avg.nc"

cntlFile=addfile(cntlFilepath,"r")
exptFile=addfile(exptFilepath,"r")

cGS=cntlFile->GROWSEAS
eGS=exptFile->GROWSEAS
dims=dimsizes(dimsizes(eGS))
if (dims.eq.3) then
  dGS=eGS(0,:,:) ; get metadata
  dGS=dGS-cGS(0,:,:)
else
  dGS=eGS ; get metadata
  dGS=dGS-cGS
end if

;************************************************
; create plot
;************************************************
;wks = gsn_open_wks("eps", "$IMG_PATH/nuclearGrowingSeason")                ; open an eps file
wks = gsn_open_wks("pdf", "$IMG_PATH/nuclearGrowingSeason")                ; open a pdf file

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
;res@cnLevels = (/-40,-30,-20,-10,-5,-2, 2,5,10,20,30,40/)
res@cnLevels = (/-60,-50,-40,-30,-20, -10, 10,20/)

res@gsnSpreadColors       = True     ; use full range of color map

res@lbLabelStride         = 1
res@lbOrientation         = "Vertical"

res@pmTickMarkDisplayMode = "Always"; use NCL default lat/lon labels

;  res@gsnAddCyclic          = False    ; data already has cyclic point
				     ; this must also be set for any zoom

; note that the gsn_csm_*map_ce templates automatically set 
; res@mpLimitMode="LatLon" for you. If you are plotting a different projection,
; you may have to set this resource.

res@mpMinLatF            = 10      ; range to zoom in on
res@mpMaxLatF            = 77

;gsn_define_colormap(wks,"BlWhRe")                  ; change colormap 
; Change color map using named colors
  colors = (/"white","black","darkorchid","blue","dodgerblue3", \
             "dodgerblue1","cyan3","cyan","white", \
             "yellow","orange"/)

  gsn_define_colormap(wks,colors)                  ; change colormap 


res@gsnLeftString        = " "
res@gsnRightString       = " "
res@gsnCenterString      = " "
res@tmXBLabelsOn         = False
res@lbLabelBarOn         = False
res@tmXBMajorLengthF     = 0.0
plot(0) = gsn_csm_contour_map_ce(wks,dGS, res)

res@mpMinLatF            = -77      ; range to zoom in on
res@mpMaxLatF            = -10
res@tmXBLabelsOn         = True
res@lbLabelBarOn         = True
res@lbOrientation        = "Horizontal"
res@lbTopMarginF         = -1.9
res@lbBottomMarginF      = 1.7
res@lbLeftMarginF        = 4.0
res@lbRightMarginF       = -1.0
res@lbLabelFontHeightF   = 0.01
;res@lbPerimOn            = True
;delete(res@tmXBMajorLengthF)
plot(1) = gsn_csm_contour_map_ce(wks,dGS, res)

pres                        = True
pres@gsnMaximize            = True
pres@gsnPaperOrientation    = "landscape"
gsn_panel(wks, plot, (/ 2, 1 /), pres)

