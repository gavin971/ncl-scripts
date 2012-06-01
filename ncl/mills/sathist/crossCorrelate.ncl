load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("plotZremaxvZmax")
procedure plotZremaxvZmax(CAM_RUN:string, summaryFile:string)

begin

  print("")
  print("k="+k+", reading data...")

  ; Get the osiris and WACCM data
  waccmRun="f40.Nabro.2deg.wcm.sd.carma.sulf.005"
  rootdir="/Volumes/Data/Models/cesm/run/"+waccmRun+"/"
  
  osirisFile = addfile(rootdir+"Osiris/osiris_hs_nabro_20110601_20110822.nc","r")
  extOsiris  = osirisFile->AEXT
  exterrOsiris  = osirisFile->AEXTerr
  lev = osirisFile->lev
  date = osirisFile->date
  preNabro=ind(date.lt.20110613)
  
  waccmFile = addfile(rootdir+"hs/f40.Nabro.2deg.wcm.sd.carma.sulf.005.cam.hs.OsirisLevs.nc","r")
  extWaccm  = waccmFile->PURSULEX750
  
  dims=dimsizes(extOsiris)
  ncol=dims(0)
  nlev=dims(1)  

  xdata = new((/ 2, ncol /), "float")
  ydata = new((/ 2, ncol /), "float")

  xdata(0, :) = sofie->Zmax(is)
;  xdata(1, :) = waccm->Zmax(iw)
  xdata(1, :) = waccm->Z532max(iw)
  
;  ydata(0, :) = sofie->NDzmax(is)
  ydata(1, :) = waccm->Zremax(iw)
  
  xdata = where(xdata.ne.0., xdata, xdata@_FillValue)
  ydata = where(ydata.ne.0., ydata, ydata@_FillValue)
  
  ; Color by particle size.
  re = wRe(iw)
  s1 = ind(re .lt. 30.)
  s2 = ind(re .ge. 30. .and. re .lt. 40.)
  s3 = ind(re .ge. 40. .and. re .lt. 50.)
  s4 = ind(re .ge. 50. .and. re .lt. 60.)
  s5 = ind(re .ge. 60.)
  

  ; Plot the correlation of the data.
  ;
  ; NOTE: Since the weather is different, this may not be very good.
  wks  = gsn_open_wks("pdf", "images/" + CAM_RUN + "/sofie_ccr_zremax_z532max")    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow")  ; choose a colormap

  res                   = True
  res@gsnDraw           = False
  res@gsnFrame          = False
;  res@tiMainString      = "Altitude Max Radius vs. Altitude Max Extinction"
  res@xyMarkLineModes   = (/ "Markers", "Markers" /)                ; choose which have markers
  res@xyMarkers         = (/ 16, 16 /)                      ; choose type of marker  
;  res@xyMarkerColors    = (/ "red", "blue" /)     ; Marker color
  res@xyMarkerSizes     = (/ 0.0075, 0.0075 /)                     ; Marker size (default 0.01)
  res@tiXAxisString     = "Zmax (km)"
  res@tiYAxisString     = "Zremax (km)"
;  res@trYLog            = True
  res@trXMinF           = 78
  res@trXMaxF           = 92
  res@trYMinF           = 78
  res@trYMaxF           = 92

  res2                  = True
  res2@gsnDraw          = False
  res2@gsnFrame         = False
  
  res@xyMarkerColors    = (/ "blue" /)     ; Marker color
  plot = gsn_csm_xy(wks, xdata(1, s1), ydata(1, s1), res)       ; contour the variable

  res@xyMarkerColors    = (/ "cyan" /)     ; Marker color
  overlay(plot, gsn_csm_xy(wks, xdata(1, s2), ydata(1, s2), res))
  
  res@xyMarkerColors    = (/ "black" /)     ; Marker color
  overlay(plot, gsn_csm_xy(wks, xdata(1,s3), ydata(1,s3), res))

  res@xyMarkerColors    = (/ "orange" /)     ; Marker color
  overlay(plot, gsn_csm_xy(wks, xdata(1,s4), ydata(1,s4), res))

  res@xyMarkerColors    = (/ "red" /)     ; Marker color
  overlay(plot, gsn_csm_xy(wks, xdata(1,s5), ydata(1,s5), res))

  overlay(plot, gsn_csm_xy(wks, (/ 0, 100. /), (/ 0, 100. /), res2))
  
    legend = create "Legend" legendClass wks 
    "vpXF"                     : 0.2                    ; orientation on page
    "vpYF"                     : 0.775
;    "vpWidthF"                 : 0.075                     ; width
;    "vpHeightF"                : 0.05                     ; height
    "vpWidthF"                 : 0.125                     ; width
    "vpHeightF"                : 0.125                     ; height
    "lgPerimOn"                : False                   ; no perimeter
    "lgItemType"               : "Markers"
    "lgMonoMarkerIndex"        : True
    "lgMarkerIndex"            : 16
    "lgItemCount"              : 5                      ; how many
    "lgLabelStrings"           : (/ " r >= 60", " 50 <= r < 60", " 40 <= r < 50", " 30 <= r < 40", " r < 30" /)
    "lgLabelsOn"               : True                   ; no default lables
;    "lgMonoLineThickness"      : False
;    "lgLineThicknesses"        : (/ 4.0, 4.0, 4.0, 4.0, 4.0 /)             ; line paters
;    "lgDashIndexes"            : (/ 0, 0, 0, 0, 0 /)
    "lgMarkerColors"             : (/ "red", "orange", "black", "cyan", "blue" /)
;    "lgMonoLineLabelFontColor" : True                    ; one label color
    "lgAutoManage"             : False
;    "lgLabelPosition"          : "Right"
    "lgLabelJust"              : "CenterLeft"
;    "lgLabelFontHeightF"       : 0.0125
    "lgLabelFontHeightF"       : 0.015
  end create

  draw(legend)  

  draw(wks)
  frame(wks)
end

