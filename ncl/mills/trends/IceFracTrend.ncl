load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;case1="b40.1955-2005.2deg.wcm.002"
;rootdir1="/data3/mmills/"+case1+"/proc/"
;rootdir2=rootdir1
;case2="b40.rcp4_5.2deg.wcm.001"
;name="Ens1"

case1="b40.1955-2005.2deg.wcm.003"
rootdir1="/data5/mmills/"+case1+"/proc/"
case2="b40.rcp4_5.2deg.wcm.002"
rootdir2="/data5/mmills/"+case2+"/proc/"
name="Ens2"

a=addfile(rootdir1+case1+".cam2.h0.ICEFRAC.MAM.nc","r")
b=addfile(rootdir2+case2+".cam2.h0.ICEFRAC.MAM.nc","r")

nyears=2007-1979+1

lat=a->lat
lon=a->lon
nlat=dimsizes(lat)
nlon=dimsizes(lon)

ICEFRAC=new((/nyears, nlat, nlon/), float)
ICEFRAC(0:2004-1979,:,:)=a->ICEFRAC(1979-1960:2004-1960,:,:)
ICEFRAC(2005-1979:2007-1979,:,:)=b->ICEFRAC(0:2,:,:)
ICEFRAC&time(26)=ICEFRAC&time(25)+365
ICEFRAC&time(27)=ICEFRAC&time(26)+365
ICEFRAC&time(28)=ICEFRAC&time(27)+365

time=ICEFRAC&time
year=time/365.d0

; regCoef calculates the linear regression coefficient between two variables.
  rc = regCoef(year,ICEFRAC(lat|:,lon|:,time|:))
  rc!0=ICEFRAC!1
  rc!1=ICEFRAC!2
  rc&lat=lat
  rc&lon=lon

; Convert from fraction to %  
  rc=rc*100.d0
  
  pdffile="~/images/IceFracTrend_"+name
  print("Creating PDF file "+pdffile+".pdf")
  wks=  gsn_open_wks ("pdf", pdffile)
;  gsn_define_colormap(wks, "wgne15")
;  gsn_define_colormap(wks, "amwg")
  colors = (/"white","black",\
             "royalblue4","royalblue3","royalblue2","royalblue1","sky blue","powder blue",\
             "dark green","green4","green3","green2","green1",\
             "seagreen3","seagreen2","seagreen1","pale green","greenyellow",\
             "lightyellow","yellow1","yellow4",\
	     "lightpink","lightpink1","lightpink2","lightpink3","lightpink4",\
             "salmon","red1","red2","red3","red4"/)

  gsn_define_colormap(wks, colors)  
  
  res                      = True               ; plot mods desired
  res@tiMainString         = "Sea Ice Extent Trends (%/year), MAM 1979-2007, "+name    ; main title
  res@tiMainFontHeightF    = 0.015
  res@cnFillOn             = True               ; turn on color fill
  res@gsnMaximize          = True          ; expand plot size to maximum
  res@lbLabelAngleF        = 60
  
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnMinLevelValF=-1.5
  res@cnMaxLevelValF=0.3
  res@cnLevelSpacingF     = .1
  res@gsnPolar   = "SH"                          ; specify the hemisphere
  res@mpMaxLatF            = -50                 ; specify max lat

; Create and draw a contour plot over a polar stereographic map.
  plot = gsn_csm_contour_map_polar(wks,rc,res)   ; plot model data  
  
  print("gv "+pdffile+".pdf")
  print("mutt -a "+pdffile+".pdf mikejmills@mac.com < /dev/null")

