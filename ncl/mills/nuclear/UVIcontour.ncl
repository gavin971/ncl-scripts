load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

path="$MODEL_PATH/ccsm/case/"

cntlCase="b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678"
exptCase="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"

cntlFilepath=path+cntlCase+"/uvi/"+cntlCase+".O3Col.h0.nc"
exptFilepath=path+exptCase+"/uvi/"+exptCase+".O3Col.h0.nc"

cntlFile=addfile(cntlFilepath,"r")
exptFile=addfile(exptFilepath,"r")

cntlUVI=cntlFile->uvi
exptUVI=exptFile->UVI
printVarSummary(exptUVI)

cntlDate=cntlFile->date
exptDate=exptFile->date

dateList=(/"20150701","20160101"/)
monthList=(/"June","December"/)

minLatList=(/-30.,-80./)
maxLatList=(/80.,30./)
;************************************************
; create plot
;************************************************
;wks = gsn_open_wks("eps", "$IMG_PATH/UVIcontour")                ; open a ps file
wks = gsn_open_wks("pdf", "$IMG_PATH/UVIcontour")                ; open a ps file

nPanel=3
plot = new(6, graphic)

do i = 0, dimsizes(dateList)-1
;i=0
  ct=ind(cntlDate.eq.dateList(i))
  et=ind(exptDate.eq.dateList(i))
  res                       = True     ; plot mods desired
  res@gsnDraw   = False
  res@gsnFrame  = False
  res@vpWidthF  = 1.0		; set width of plot 
  res@vpHeightF = 1.0		; set height of plot

  res@tmXBMajorLengthF     = 0.0
  res@cnFillOn              = True     ; turn on color fill
  res@cnLinesOn             = True    ; turn of contour lines
;  res@cnLevelSelectionMode  = "ExplicitLevels"
;  res@cnLevels = (/2.5,5.5,7.5,10.5,13.5,16.5,19.5,22.5,25.5/)
  res@cnLevelSelectionMode  = "ManualLevels"
  res@cnLevelSpacingF = 3
  res@cnMinLevelValF  = 3
  res@cnMaxLevelValF  = 27

  res@gsnSpreadColors       = False     ; use full range of color map

  res@lbLabelStride         = 1
  res@lbOrientation         = "Vertical"

  res@pmTickMarkDisplayMode = "Always"; use NCL default lat/lon labels

;  res@gsnAddCyclic          = False    ; data already has cyclic point
				       ; this must also be set for any zoom

; note that the gsn_csm_*map_ce templates automatically set 
; res@mpLimitMode="LatLon" for you. If you are plotting a different projection,
; you may have to set this resource.

  res@mpMinLatF            = minLatList(i)      ; range to zoom in on
  res@mpMaxLatF            =  maxLatList(i)
;  res@mpMinLonF            =  30.
;  res@mpMaxLonF            = 120.

; Change color map using named colors
;  colors = (/"white","black", \
;             "green4","yellow", \
;             "orange","red", \
;             "purple","tan1", \
;             "tan3","tan4","black","white"/)
  colors = (/"white","black","blue","green","yellow", \
             "orange","red","red4", \
             "sienna4","purple","black","white"/)

  gsn_define_colormap(wks,colors)                  ; change colormap 

  res@lbLabelBarOn=False
  if (i.eq.0) then 
    res@lbLabelBarOn=True
  end if
;  res@lbTopMarginF=0.4
  res@lbBottomMarginF=-1.7
;  res@lbLeftMarginF=1
;  res@lbRightMarginF=50
  res@lbLabelFontHeightF   = 0.018
  res@tmXBLabelFontHeightF   = 0.018
  res@tmYLLabelFontHeightF   = 0.018
  res@gsnLeftString=" "
  res@gsnRightString=" "
;  res@gsnCenterString="UVI "+monthList(i)+", Control"
  res@gsnCenterString=" "
  res@gsnCenterStringFontHeightF = 0.02
  plot(i) = gsn_csm_contour_map_ce(wks,cntlUVI(ct,:,:), res)
  
;  res@gsnCenterString="UVI "+monthList(i)+", Year 3"
  res@gsnCenterString=" "
  res@tmXBLabelsOn         = True
  res@lbLabelBarOn=False
  plot(i+2) = gsn_csm_contour_map_ce(wks,exptUVI(et,:,:), res)
  
;  gsn_define_colormap(wks,"WhiteYellowOrangeRed")                  ; change colormap 
;  res@gsnSpreadColors       = True     ; use full range of color map
  diffUVI=exptUVI(et,:,:)
  diffUVI=diffUVI-cntlUVI(ct,:,:)
  delete(res@cnLevels)
;  res@cnLevelSelectionMode  = "ExplicitLevels"
;  res@cnLevels = (/-0.5,1,2,3,4,4.5,5,5.5,6/)
  res@cnLevelSelectionMode  = "ManualLevels"
  res@cnLevelSpacingF = 1
  res@cnMinLevelValF  = -1
  res@cnMaxLevelValF  = 6
;  res@gsnCenterString="UVI "+monthList(i)+", Year 3 - Control"
  res@gsnCenterString=" "
  res@lbBottomMarginF= 0.0
  res@tmXBLabelsOn         = False
  if (i.eq.0) then 
    res@lbLabelBarOn=True
  end if
  plot(i+4) = gsn_csm_contour_map_ce(wks,diffUVI(:,:), res)
  delete(res@cnLevels)
    
end do

;lbres                           = True
;lbres@lbAutoManage              = False         ; we control label bar
;lbres@lbFillColors              = fill_colors   ; use nice strided colors  
;lbres@lbPerimOn                 = False         ; Turn off labelbar perimeter.
;lbres@lbMonoFillPattern         = True          ; one pattern, all solid
;lblabels    = "" + fspan(3,27,9)

;nboxes=8
;gsn_labelbar_ndc(wks,nboxes,lblabels,0.5,0.2,lbres)

pres                        = True
pres@gsnMaximize            = True
pres@gsnPaperOrientation    = "landscape"
gsn_panel(wks, plot, (/ 3, 2 /), pres)
;frame(wks)
