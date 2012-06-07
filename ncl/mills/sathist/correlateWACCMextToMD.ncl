; Comparison of WACCM Nabro extinction calculations with Osiris observations

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("correlateWACCMtoOsiris")
procedure correlateWACCMtoOsiris(k:integer)

begin

  print("")
  print("k="+k+", reading data...")

  ; Get the osiris and WACCM data
  waccmRun="f40.Nabro.2deg.wcm.sd.carma.sulf.005"
  rootdir="/Volumes/Data/Models/cesm/run/"+waccmRun+"/"

  waccmFile = addfile(rootdir+"hs/f40.Nabro.2deg.wcm.sd.carma.sulf.005.cam.hs.nc","r")
  ext  = waccmFile->PURSULEX750
  sulfmd  = waccmFile->PURSULMD
  z3  = waccmFile->Z3
  zmin=min(z3(:,k))/1000.0
  zmax=max(z3(:,k))/1000.0
    
  lev = waccmFile->lev
  date = waccmFile->date
  preNabro=ind(date.lt.20110613)
  week1=ind(date.ge.20110613.and.date.lt.20110620)
  week2=ind(date.ge.20110620.and.date.lt.20110627)
  week3=ind(date.ge.20110627.and.date.lt.20110704)
  week4=ind(date.ge.20110704.and.date.lt.20110711)
  postweek4=ind(date.ge.20110711)
  
  dims=dimsizes(ext)
  ncol=dims(0)
  nlev=dims(1)
  
  data = new((/ 2, ncol /), "float")
  data(0, :) = ext(:,k)
  data(1, :) = (/sulfmd(:,k)/)
  data@_FillValue = -999.
  
  ; Plot the correlation of the data.
  pdfpath="$IMG_PATH/"+waccmRun+"/Osiris/"
  a = systemfunc("mkdir -p "+pdfpath)
  pdffile=pdfpath+"corr_ext_sulfmd_"+sprintf("%4.1f",lev(k))+"hPa"
  print("Creating file "+pdffile+".pdf")
  wks  = gsn_open_wks("pdf", pdffile)    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow+white+gray")  ; choose a colormap

  res                   = True
  res@gsnDraw           = False
  res@gsnFrame          = False
  res@tiMainString      = "WACCM/CARMA at "+sprintf("%4.1f",lev(k))+" hPa, "+\
    sprintf("%4.1f",zmin)+"-"+sprintf("%4.1f",zmax)+" km"
  res@xyMarkLineModes   = "Markers"                ; choose which have markers
  res@xyMarkers         =  16                      ; choose type of marker  
  res@xyMarkerColor     = "red"                    ; Marker color
  res@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
  res@tiXAxisString     = "Extinction (750 nm) (/km)"
  res@tiYAxisString     = "Sulfate Mass Density (g/cm3)"
  res@trXMinF           = 0.0
  res@trYMinF           = 0.0
;  res@trXMaxF           = max(data)
;  res@trYMaxF           = max(data)
  
  print("Plotting post week 4 scatterplot...")
  plot = gsn_csm_xy(wks, ext(postweek4,k), sulfmd(postweek4,k), res)       ; scatterplot
  
  print("Plotting week 4 scatterplot...")
  res@xyMarkerColor     = "orange"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, ext(week4,k), sulfmd(week4,k), res))       ; scatterplot
  
  print("Plotting week 3 scatterplot...")
  res@xyMarkerColor     = "yellow"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, ext(week3,k), sulfmd(week3,k), res))       ; scatterplot
  
  print("Plotting week 2 scatterplot...")
  res@xyMarkerColor     = "green"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, ext(week2,k), sulfmd(week2,k), res))       ; scatterplot
  
  print("Plotting week 1 scatterplot...")
  res@xyMarkerColor     = "blue"                    ; Marker color
  res@xyMarkerSizeF     = 0.005                     ; Marker size (default 0.01)
  overlay(plot, gsn_csm_xy(wks, ext(week1,k), sulfmd(week1,k), res))       ; scatterplot
  
  print("Plotting pre-Nabro scatterplot...")
  res@xyMarkerColor     = "gray"                    ; Marker color
  res@xyMarkerSizeF     = 0.002                     ; Marker size (default 0.01)
  overlay(plot, gsn_csm_xy(wks, ext(preNabro,k), sulfmd(preNabro,k), res))       ; scatterplot
    
  draw(plot)
  frame(wks)
  
;  print("Plotting geopotential height...")
;  res@gsnDraw           = True
;  res@gsnFrame          = True
;  res@tiMainString      = "Geopotential height at "+sprintf("%4.1f",lev(k))+" hPa"
;  res@xyMarkerColor     = "black"                    ; Marker color
;  delete(res@trYMinF)
;  res@tiYAxisString     = "Geopotential height (km)"
;  
;  plot = gsn_csm_xy(wks, ext(:,k), z3(:,k)/1000., res)       ; scatterplot
  
  delete(wks)
  delete(plot)
  delete(res)

end

;correlateWACCMtoOsiris(54) ; 118 hPa
;correlateWACCMtoOsiris(53) ; 100.5 hPa
;correlateWACCMtoOsiris(52) ; 85.4 hPa
;correlateWACCMtoOsiris(51) ; 72.6 hPa
;correlateWACCMtoOsiris(50) ; 61.5 hPa
;correlateWACCMtoOsiris(49) ; 52.0 hPa
correlateWACCMtoOsiris(48) ; 43.9 hPa
