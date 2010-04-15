load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

sootrun="vgsoot08"
controlrun="vgcontrol02"

rootdir="/data5/mmills/"

soot=addfile(rootdir+sootrun+"/concat_tem_diag0_"+sootrun+"_y2045-2064_tropics_avg.nc","r")
control=addfile(rootdir+controlrun+"/concat_tem_diag0_"+controlrun+"_y2045-2064_tropics_avg.nc","r")

wstar_soot=soot->WSTAR
wstar_cntl=control->WSTAR
lev=control->lev

data = new((/2,dimsizes(wstar_soot)/),double)

data(0,:)=wstar_soot*1000.d0 ; convert m/s to mm/s
data(1,:)=wstar_cntl*1000.d0


;************************************************
; plotting parameters
;************************************************
  wks   = gsn_open_wks ("pdf","~/images/vgsoot_wstar_tropics")      ; open workstation
  res                   = True                     ; plot mods desired
  res@gsnPaperOrientation = "portrait"
  res@gsnMaximize       = True
  res@gsnBoxMargin      = 0.0
  res@xyDashPattern = 0               ; Make curves all solid
  res@xyLineColors = (/"Red","Blue"/)
  
  res@trYLog            = True                       ; log Y-axis
  res@trYMinF           = 3
  res@trYMaxF           = 600.0
  res@trYReverse        = True                       ; reverse Y-axis

  res@trXMinF           = 0.0
  res@trXMaxF           = 1.5
;  res@trXLog            = True                       ; log X-axis
  
  res@tiYAxisString     ="Pressure (hPa)"
  res@tiXAxisString     ="Vertical winds (mm/s)"
  res@tiMainString      = "Tropics Average"             ; add title
  
; Add altitude tick labels on the right Y-axis:  
  add_hgt=True
  set_pres_hgt_axes(lev,res,add_hgt)
 
  plot  = gsn_csm_xy (wks,data,lev,res) ; create plot

;
; Add a right Y axis "height" label.
;
  if(add_hgt)
    rightaxis_string = create "right_axis" textItemClass wks
      "txString"      : "Height (km)"
      "txFontHeightF" : 0.025
      "txAngleF"      : 90.
    end create

    anno = NhlAddAnnotation(plot,rightaxis_string)

    setvalues anno
      "amZone"          : 3      ; Just outside plot area
      "amJust"          : "centercenter"
      "amSide"          : "right" 
      "amParallelPosF"  : 0.5
      "amOrthogonalPosF": 0.03
      "amResizeNotify"  : True     ; Resize if plot resized.
    end setvalues
  end if
