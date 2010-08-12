load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;  run1="b40.1850.track1.2deg.wcm.007"
  run1="b40.1850.2deg.wset.002"
  tspath="/data3/mmills/"+run1+"/h0.annual/"

  historyPattern=tspath+"*.TS.annual.nc"
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  print("historyNames:"+dimsizes(historyNames))

  lat = history[0]->lat
  lon = history[0]->lon
  print("Reading time...")
  time = addfiles_GetVar(history, historyNames, "time")
  print("Reading surface temperature...")
  ts = addfiles_GetVar(history, historyNames, "TS")
  
  decade = time/3650.d0
  decade@units="decade"
  rc = regCoef(decade,ts(lat|:,lon|:,time|:))
  rc!0=ts!1
  rc!1=ts!2
  rc&lat=lat
  rc&lon=lon
  
  pdffile="~/images/"+run1+"_TStrend"
  print("Creating PDF file "+pdffile+".pdf")
  wks=  gsn_open_wks ("pdf", pdffile)
  gsn_define_colormap(wks, "wgne15")
  
  res                      = True               ; plot mods desired
  res@tiMainString         = "Surface Temperature Trends (K/decade), "+run1    ; main title
  res@tiMainFontHeightF    = 0.015
  res@cnFillOn             = True               ; turn on color fill
  res@gsnMaximize          = True          ; expand plot size to maximum
  res@lbLabelAngleF        = 60
  
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnMinLevelValF=-.36
  res@cnMaxLevelValF=.28
  res@cnLevelSpacingF     = .04
  
  plot = gsn_csm_contour_map_ce(wks,rc,res)   ; plot model data  
  
  print("gv "+pdffile+".pdf")
  print("mutt -a "+pdffile+".pdf mikejmills@mac.com < /dev/null")

