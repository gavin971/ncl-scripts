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
wks = gsn_open_wks("eps", "$IMG_PATH/nuclearGrowingSeason")                ; open an eps file
;wks = gsn_open_wks("pdf", "$IMG_PATH/nuclearGrowingSeason")                ; open a pdf file

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

res@mpMinLatF            = 15      ; range to zoom in on
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

res@mpMinLatF            = -55      ; range to zoom in on
res@mpMaxLatF            = -15
res@tmXBLabelsOn         = True
res@lbLabelBarOn         = True
res@lbOrientation        = "Horizontal"
res@lbTopMarginF         = 1.0
res@lbBottomMarginF      = -1.2
res@lbLabelFontHeightF   = 0.01
delete(res@tmXBMajorLengthF)
plot(1) = gsn_csm_contour_map_ce(wks,dGS, res)

;getvalues plot(1)                      ; Get plot size for use in
;  "vpHeightF" : vph                 ; creating labelbar.
;  "vpWidthF"  : vpw
;end getvalues

;nboxes = dimsizes(colors)

;lbres                    = True          ; labelbar only resources
;lbres@lbAutoManage       = True          ; Necessary to control sizes
;lbres@vpWidthF           = 0.4 * vpw     ; labelbar width
;lbres@vpHeightF          = 0.5 * vph     ; labelbar height
;lbres@lbBoxMajorExtentF  = 0.0          ; puts space between color boxes
;lbres@lbFillColors       = colors ; labelbar colors
;lbres@lbMonoFillPattern  = True          ; Solid fill pattern
;lbres@lbLabelFontHeightF = 0.08          ; font height. default is small
;lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
;lbres@lbOrientation      = "Horizontal"
;lbres@lbPerimOn          = False

;labels = res@cnLevels

;lbid = gsn_create_labelbar(wks,nboxes,labels,lbres)
   
;
; Now, create some annotation resources indicating how we want to
; attach the labelbar to the plot. Here, we are using the top right
; corner of the labelbar as the point which we are going to position
; it, and then we use amParallelPosF and amOrthogonalPosF to indicate
; where we want to place it.
;
; amParallelPosF/amOrthogonalPosF                                       ;
;    0.0/ 0.0  -  annotation in dead center of plot
;    0.5/ 0.5  -  annotation at bottom right of plot
;    0.5/-0.5  -  annotation at top right of plot
;   -0.5/-0.5  -  annotation at top left of plot
;   -0.5/ 0.5  -  annotation at bottom left of plot
;                                                                 
;amres                  = True
;amres@amJust           = "TopRight"
;amres@amParallelPosF   =  0.5
;amres@amOrthogonalPosF =  0.2
;annoid = gsn_add_annotation(plot(1),lbid,amres)

pres                        = True
pres@gsnMaximize            = True
pres@gsnPaperOrientation    = "landscape"
gsn_panel(wks, plot, (/ 2, 1 /), pres)


