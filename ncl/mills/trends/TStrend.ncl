load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

  host=getenv("HOSTNAME")
  run1="b40.1955-2005.2deg.wcm.avg"
;  run1="b40.1955-2005.2deg.wset.avg"

  if (host.eq."cmip5") then
    if (isStrSubset(run1,"wcm")) then
      tspath = "/data/waccm/"
    else
      tspath = "/data/ccsm4/"
    end if
  else
    tspath = "/data3/mmills/"
  end if

  tspath=tspath+run1+"/h0.annual/"

;  year1=1986
;  year2=2005
  
;  year1=1960
;  year2=1979

  year1=1960
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
  
  pdffile="$IMAGE_PATH/"+run1+"_TStrend_"+year1+"-"+year2
  print("Creating PDF file "+pdffile+".pdf")
  wks=  gsn_open_wks ("pdf", pdffile)
  gsn_define_colormap(wks, "wgne15")
  
  res                      = True               ; plot mods desired
  res@tiMainString         = "Surface Temperature Trends (K/decade), "+run1+":"+year1+"-"+year2    ; main title
  res@tiMainFontHeightF    = 0.015
  res@cnFillOn             = True               ; turn on color fill
  res@gsnMaximize          = True          ; expand plot size to maximum
  res@lbLabelAngleF        = 60
  
  res@cnLevelSelectionMode = "ManualLevels" 
;  res@cnMinLevelValF=-2.8
;  res@cnMaxLevelValF=2.8
;  res@cnLevelSpacingF     = .4
  res@cnMinLevelValF=-2.1
  res@cnMaxLevelValF=2.1
  res@cnLevelSpacingF     = .3

; Create and draw a contour plot over a cylindrical equidistant map.  
  plot = gsn_csm_contour_map_ce(wks,rc,res)   ; plot model data  
  
  print("gv "+pdffile+".pdf")
  print("mutt -a "+pdffile+".pdf -s "+pdffile+" mikejmills@mac.com < /dev/null")

