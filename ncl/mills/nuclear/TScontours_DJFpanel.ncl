; Script written by Mike Mills
; Modified September 24, 2013

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;path="$MODEL_PATH/cesm/runs/cesm1_0/"
path="/glade/p/acd/mmills/run/"

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
;wks = gsn_open_wks("eps", "$IMG_PATH/nuclearTSdiff")    ; open an eps file
;wks = gsn_open_wks("pdf", "$IMG_PATH/nuclear/Fig10_TSdiff_hotcold")                ; open a pdf file
wks = gsn_open_wks("png", "$IMG_PATH/nuclear/TSdiff_DJF")                ; open a pdf file

;nPanel=2
nPanel=1
plot = new(nPanel, graphic)

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
;  colors = (/"white","black","darkorchid","orchid","blue", \
;             "dodgerblue1","cyan4","cyan1","seagreen1","green4", \
;             "green","greenyellow","white", "lightgoldenrod","orange"/)

;  colors = (/"blue4","blue","aquamarine","white","orange"/)

;  rgb_array = span_named_colors(colors,False)

;  gsn_define_colormap(wks,rgb_array)
;  gsn_draw_colormap(wks)

;---Using the color map to define a color map for the workstation
; gsn_define_colormap(wks,"testcmap")
 gsn_define_colormap(wks,"hotcold_18lev")
; gsn_define_colormap(wks,"BlueDarkRed18") 
 
res@gsnSpreadColors = True
;res@gsnSpreadColorStart = 0
res@gsnSpreadColorEnd = 13
res@gsnLeftString        = " "
res@gsnRightString       = " "
res@gsnCenterString      = " "
res@tmXBLabelsOn         = True
res@lbLabelBarOn         = True
res@lbTopMarginF         = 0.0
res@lbBottomMarginF      = -2.8
;res@lbRightMarginF         = 0.3
res@tmXBMajorLengthF     = 0.0
;plot(0) = gsn_csm_contour_map_ce(wks,dTSJJA(0,:,:), res)

txres = True
txres@txFontColor="black"
txres@txPerimOn=True
txres@txBackgroundFillColor="lightgoldenrod1"
txres@txFontHeightF = 0.03                     ; font smaller. default big
xpos=-160
ypos=0
;dum=gsn_add_text(wks,plot(0),"a",xpos,ypos,txres) 

;res@tmXBLabelsOn         = False
res@tmXBLabelsOn         = True
;res@tmXTLabelsOn         = True
;res@tmYRLabelsOn         = True
;res@tmXTOn         = True
;res@tmYROn         = True
res@tmXTMajorLengthF     = 0.0
res@tmYRMajorLengthF     = 0.0

res@lbLabelBarOn         = True
;res@lbOrientation        = "Vertical"
res@lbOrientation        = "Horizontal"
;res@lbTopMarginF         = 5
;res@lbBottomMarginF      = -7
res@lbTopMarginF         = -0.1
res@lbBottomMarginF      = 0.1
res@lbLeftMarginF         = -0.7
res@lbRightMarginF      = -0.7
res@lbLabelFontHeightF   = 0.015

;delete(res@tmXBMajorLengthF)
;plot(1) = gsn_csm_contour_map_ce(wks,dTSDJF(0,:,:), res)
plot(0) = gsn_csm_contour_map_ce(wks,dTSDJF(0,:,:), res)
;dum=gsn_add_text(wks,plot(1),"b",xpos,ypos,txres) 

;draw(plot)
;frame(wks)
;delete(wks)

pres                        = True
pres@gsnMaximize            = True
pres@gsnPaperOrientation    = "landscape"
;gsn_panel(wks, plot, (/ 2, 1 /), pres)
gsn_panel(wks, plot, (/ 1, 1 /), pres)


