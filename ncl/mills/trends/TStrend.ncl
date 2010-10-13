load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;  run1="b40.1850.track1.2deg.wcm.007"
;  run1="b40.1850.2deg.wset.002"
  run1="b40.20th.2deg.wset.001"
;  run1="b40.1955-2005.2deg.wcm.002"
  tspath="/data3/mmills/"+run1+"/h0.annual/"

  year1=1986
  year2=2005

  nyears=year2-year1+1
  historyNames=new(nyears,string)

  do year = year1, year2
    date=year*10000+701
    if ( year .lt. 10 ) then
        yearstr="000"+year
    else
      if ( year .lt. 100 ) then
        yearstr="00"+year
      else
        if ( year .lt. 1000 ) then
          yearstr="0"+year
        else
          yearstr=""+year
        end if
      end if
    end if

    historyNames(year-year1)=tspath+run1+".cam2.h0."+yearstr+".annual.nc"
    
  end do
  
  print(historyNames)

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
  
  print("decade:"+dimsizes(decade))
  print("ts:"+dimsizes(ts))
  
; regCoef calculates the linear regression coefficient between two variables.
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
  res@cnMinLevelValF=-2.8
  res@cnMaxLevelValF=2.8
  res@cnLevelSpacingF     = .4

; Create and draw a contour plot over a cylindrical equidistant map.  
  plot = gsn_csm_contour_map_ce(wks,rc,res)   ; plot model data  
  
  print("gv "+pdffile+".pdf")
  print("mutt -a "+pdffile+".pdf mikejmills@mac.com < /dev/null")

