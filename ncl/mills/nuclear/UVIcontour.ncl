load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;path="$MODEL_PATH/ccsm/case/"
;path="$MODEL_PATH/cesm/runs/cesm1_0/"
path="$MODEL_PATH/"

cntlCase="b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678"
exptCase="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"

cntlFilepath=path+cntlCase+"/uvi/"+cntlCase+".O3Col.h0.nc"
exptFilepath=path+exptCase+"/uvi/"+exptCase+".O3Col.h0.nc"

cntlFile=addfile(cntlFilepath,"r")
exptFile=addfile(exptFilepath,"r")

cntlUVI=cntlFile->uvi
exptUVI=exptFile->UVI ; upper case UVI includes adjustment for soot optical depth
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
wks_type = "eps"
wks_type@wkPaperWidthF  = 11.0  ; in inches
wks_type@wkPaperHeightF = 14.0  ; in inches
;wks = gsn_open_wks(wks_type, "$IMG_PATH/nuclear/Fig09_UVIcontour")                ; open a ps file
wks = gsn_open_wks("pdf", "$IMG_PATH/nuclear/Fig09_UVIcontour")                ; open a ps file

nPanel=3
plot = new(6, graphic)

do i = 0, dimsizes(dateList)-1
;i=0
  ct=ind(cntlDate.eq.dateList(i))
  et=ind(exptDate.eq.dateList(i))
  res                       = True     ; plot mods desired
  res@gsnDraw   = True
  res@gsnFrame  = False
;  res@vpWidthF  = 1.0		; set width of plot 
;  res@vpHeightF = 1.0		; set height of plot

  res@tmXBMajorLengthF     = 0.0
  res@cnFillOn              = True     ; turn on color fill
  res@cnLinesOn             = True    ; turn of contour lines
  res@cnLineLabelsOn        = True     ; turn on contour line labels
  res@cnLineLabelFontHeightF=0.01
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

  res@lbLabelBarOn=False
  res@tiXAxisOffsetYF = 0.005
  if (i.eq.0) then 
    res@lbLabelBarOn=True
    res@tiXAxisOffsetYF = 0.02
  end if
  res@lbOrientation        = "Horizontal"
  res@lbTopMarginF=0.5
  res@lbBottomMarginF=-0.5
  res@lbLeftMarginF=6.0
  res@lbRightMarginF=-7.0
  res@lbLabelFontHeightF   = 0.01
  res@tmXBLabelFontHeightF   = 0.01
  res@tmYLLabelFontHeightF   = 0.01
  res@gsnLeftString=" "
  res@gsnRightString=" "
;  res@gsnCenterString="UVI "+monthList(i)+", Control"
  res@gsnCenterString=" "
  res@gsnCenterStringFontHeightF = 0.02
  res@tiXAxisFontHeightF = 0.01
  res@tiXAxisString = "Longitude"
  res@tiYAxisFontHeightF = 0.01
  res@tiYAxisString = "Latitude"
  res@tmXBLabelFontHeightF = 0.01
  res@tmYRLabelFontHeightF = 0.01
  
  res@vpHeightF = 0.3 
  res@vpWidthF = 0.39 
  res@vpXF = 0.1+i*0.47
  res@vpYF = 0.99 
  ct=ind(cntlDate.eq.dateList(i))
  et=ind(exptDate.eq.dateList(i))
  
; Change color map using named colors
  colors = (/"white","black","blue","green","yellow", \
             "orange","red","red4", \
             "sienna4","purple","black","white"/)

  gsn_define_colormap(wks,colors)                  ; change colormap 

  plot(i) = gsn_csm_contour_map_ce(wks,cntlUVI(ct,:,:), res)
;  if (i.eq.0) then
;    drawNDCGrid(wks)
;  end if
  
;  res@gsnCenterString="UVI "+monthList(i)+", Year 3"
  res2=res
  res2@tiXAxisOffsetYF = 0.005
  res2@tmXBLabelsOn         = True
  res2@lbLabelBarOn=False
  res2@vpYF = 0.8 
  plot(i+2) = gsn_csm_contour_map_ce(wks,exptUVI(et,:,:), res2)
  
  colors2 = (/"white","black","lightblue","white", \
             "pink","salmon","red1","red3","red4","black","white"/)
  
  gsn_define_colormap(wks,colors2)                  ; change colormap 
  res@gsnSpreadColors       = False     ; use full range of color map
  diffUVI=exptUVI(et,:,:)
  diffUVI=diffUVI-cntlUVI(ct,:,:)
  delete(res@cnLevels)
  res3=res2
  res3@cnLevelSelectionMode  = "ExplicitLevels"
  res3@cnLevels = (/-0.5,0.5,1.5,2.5,3.5,4.5,5.5/)
;  res3@cnLevelSelectionMode  = "ManualLevels"
;  res3@cnLevelSpacingF = 1.0
;  res3@cnMinLevelValF  = -0.5
;  res3@cnMaxLevelValF  = 5.5
;  res@gsnCenterString="UVI "+monthList(i)+", Year 3 - Control"
  res3@gsnCenterString=" "
;  res3@lbBottomMarginF= 0.0
;  res3@tmXBLabelsOn         = False
  res3@tiXAxisOffsetYF = 0.005
  if (i.eq.0) then 
    res3@lbLabelBarOn=True
    res3@tiXAxisOffsetYF = 0.02
  end if
  res3@vpYF = 0.62 
  plot(i+4) = gsn_csm_contour_map_ce(wks,diffUVI(:,:), res3)
  delete(res@cnLevels)
    
end do

txres = True
txres@txFontColor="black"
txres@txPerimOn=True
txres@txBackgroundFillColor="lightgoldenrod1"
txres@txFontHeightF = 0.015                     ; font smaller. default big
xpos=-160
ypos=65
xpos2=-160
ypos2=15
gsn_text(wks,plot(0),"a",xpos,ypos,txres) 
gsn_text(wks,plot(1),"b",xpos2,ypos2,txres) 
gsn_text(wks,plot(2),"c",xpos,ypos,txres) 
gsn_text(wks,plot(3),"d",xpos2,ypos2,txres) 
gsn_text(wks,plot(4),"e",xpos,ypos,txres) 
gsn_text(wks,plot(5),"f",xpos2,ypos2,txres) 

;lbres                           = True
;lbres@lbAutoManage              = False         ; we control label bar
;lbres@lbFillColors              = fill_colors   ; use nice strided colors  
;lbres@lbPerimOn                 = False         ; Turn off labelbar perimeter.
;lbres@lbMonoFillPattern         = True          ; one pattern, all solid
;lblabels    = "" + fspan(3,27,9)

;nboxes=8
;gsn_labelbar_ndc(wks,nboxes,lblabels,0.5,0.2,lbres)

;pres                        = True
;pres@gsnMaximize            = True
;pres@gsnPaperOrientation    = "landscape"
;gsn_panel(wks, plot, (/ 3, 2 /), pres)
frame(wks)
;psres = True   
;maximize_output(wks,psres) 
