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
wks = gsn_open_wks("eps", "$IMG_PATH/UVIcontour")                ; open a ps file
;wks = gsn_open_wks("pdf", "$IMG_PATH/UVIcontour")                ; open a ps file

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

  res@cnFillOn              = True     ; turn on color fill
  res@cnLinesOn             = True    ; turn of contour lines
  res@cnLevelSelectionMode  = "ExplicitLevels"
  res@cnLevels = (/2.5,5.5,7.5,10.5,13.5,16.5,19.5,22.5,25.5/)

  res@gsnSpreadColors       = False     ; use full range of color map

  res@lbLabelStride         = 1

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
  colors = (/"white","black", \
             "green4","yellow", \
             "orange","red", \
             "purple","sienna1", \
             "cadetblue","darkorange4","midnightblue","black"/)

  gsn_define_colormap(wks,colors)                  ; change colormap 

  res@lbLabelBarOn=False
  if (i.eq.0) then 
    res@lbLabelBarOn=True
  end if
  res@lbTopMarginF=0.4
;  res@lbLeftMarginF=1
;  res@lbRightMarginF=50
  res@gsnLeftString=" "
  res@gsnRightString=" "
  res@gsnCenterString="UVI "+monthList(i)+", Control"
  plot(i) = gsn_csm_contour_map_ce(wks,cntlUVI(ct,:,:), res)
  
  res@gsnCenterString="UVI "+monthList(i)+", Year 3"
  res@lbLabelBarOn=False
  plot(i+2) = gsn_csm_contour_map_ce(wks,exptUVI(et,:,:), res)
  
;  gsn_define_colormap(wks,"WhiteYellowOrangeRed")                  ; change colormap 
;  res@gsnSpreadColors       = True     ; use full range of color map
  diffUVI=exptUVI(et,:,:)
  diffUVI=diffUVI-cntlUVI(ct,:,:)
  delete(res@cnLevels)
  res@cnLevels = (/0,1,2,3,4,4.5,5,5.5,6/)
  res@gsnCenterString="UVI "+monthList(i)+", Year 3 - Control"
  if (i.eq.0) then 
    res@lbLabelBarOn=True
  end if
  plot(i+4) = gsn_csm_contour_map_ce(wks,diffUVI(:,:), res)
  delete(res@cnLevels)
    
end do

pres                        = True
pres@gsnMaximize            = True
pres@gsnPaperOrientation    = "landscape"
gsn_panel(wks, plot, (/ 3, 2 /), pres)
