load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

path="$MODEL_PATH/ccsm/case/"

cntlCase="b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678"
exptCase="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"


cntlFilepath=path+cntlCase+"/proc/"+cntlCase+".cam.growseas.2014-2018avg.nc"
exptFilepath=path+exptCase+"/proc/"+exptCase+".cam.growseas.2014-2018avg.nc"

cntlFile=addfile(cntlFilepath,"r")
exptFile=addfile(exptFilepath,"r")

cGS=cntlFile->GROWSEAS
eGS=exptFile->GROWSEAS
dGS=eGS ; get metadata
dGS=dGS-cGS

;************************************************
; create plot
;************************************************
wks = gsn_open_wks("pdf", "$IMG_PATH/nuclearGrowingSeason")                ; open a ps file

;nPanel=2
;plot = new(2, graphic)

res                       = True     ; plot mods desired
;res@gsnDraw   = False
;res@gsnFrame  = False
;res@vpWidthF  = 1.0		; set width of plot 
;res@vpHeightF = 1.0		; set height of plot

res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = True    ; turn of contour lines
res@cnLevelSelectionMode  = "ExplicitLevels"
res@cnLevels = (/-50,-40,-30,-20,-10,-5,-2, 2,5,10,20,30,40,50/)

res@gsnSpreadColors       = True     ; use full range of color map

res@lbLabelStride         = 1
res@lbOrientation         = "Vertical"

res@pmTickMarkDisplayMode = "Always"; use NCL default lat/lon labels

;  res@gsnAddCyclic          = False    ; data already has cyclic point
				     ; this must also be set for any zoom

; note that the gsn_csm_*map_ce templates automatically set 
; res@mpLimitMode="LatLon" for you. If you are plotting a different projection,
; you may have to set this resource.

res@mpMinLatF            = -60      ; range to zoom in on
res@mpMaxLatF            = 75

gsn_define_colormap(wks,"BlWhRe")                  ; change colormap 

res@gsnLeftString=" "
res@gsnRightString=" "
res@gsnCenterString=" "
plot = gsn_csm_contour_map_ce(wks,dGS, res)

;res@mpMinLatF            = -80      ; range to zoom in on
;res@mpMaxLatF            = -10

;plot(1) = gsn_csm_contour_map_ce(wks,dGS, res)


;pres                        = True
;pres@gsnMaximize            = True
;pres@gsnPaperOrientation    = "landscape"
;gsn_panel(wks, plot, (/ 2, 1 /), pres)


