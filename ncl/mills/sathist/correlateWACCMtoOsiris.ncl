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
  
  osirisFile = addfile(rootdir+"Osiris/osiris_hs_nabro_20110601_20110822.nc","r")
  extOsiris  = osirisFile->AEXT
  exterrOsiris  = osirisFile->AEXTerr
  lev = osirisFile->lev
  date = osirisFile->date
  preNabro=ind(date.lt.20110613)
  week1=ind(date.ge.20110613.and.date.lt.20110620)
  week2=ind(date.ge.20110620.and.date.lt.20110627)
  week3=ind(date.ge.20110627.and.date.lt.20110704)
  week4=ind(date.ge.20110704.and.date.lt.20110711)
  postweek4=ind(date.ge.20110711)
  
  waccmFile = addfile(rootdir+"hs/f40.Nabro.2deg.wcm.sd.carma.sulf.005.cam.hs.OsirisLevs.nc","r")
  extWaccm  = waccmFile->PURSULEX750
  
  extOsiris  = where(ismissing(extOsiris),0.,extOsiris)
  
  dims=dimsizes(extOsiris)
  ncol=dims(0)
  nlev=dims(1)
  
  data = new((/ 2, ncol /), "float")
  data(0, :) = extWaccm(:,k)
  data(1, :) = (/extOsiris(:,k)/)
  data@_FillValue = -999.
  
  ; Plot the correlation of the data.
  pdfpath="$IMG_PATH/"+waccmRun+"/Osiris/"
  a = systemfunc("mkdir -p "+pdfpath)
  pdffile=pdfpath+"corr_ext_"+sprintf("%04.1f",lev(k)/1000.)+"km"
  print("Creating file "+pdffile+".pdf")
  wks  = gsn_open_wks("pdf", pdffile)    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow+white+gray")  ; choose a colormap

  res                   = True
  res@gsnDraw           = False
  res@gsnFrame          = False
  res@tiMainString      = "Extinction (750 nm) (/km) at "+sprintf("%4.1f",lev(k)/1000.)+" km"
  res@xyMarkLineModes   = "Markers"                ; choose which have markers
  res@xyMarkers         =  16                      ; choose type of marker  
  res@xyMarkerColor     = "red"                    ; Marker color
  res@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
  res@tiXAxisString     = "WACCM/CARMA model"
  res@tiYAxisString     = "Osiris observations"
  res@trXMinF           = 0.0
  res@trYMinF           = 0.0
  res@trXMaxF           = max(data)
  res@trYMaxF           = max(data)
  
  print("Plotting post week 4 scatterplot...")
  plot = gsn_csm_xy(wks, extWaccm(postweek4,k), extOsiris(postweek4,k), res)       ; scatterplot
  
  print("Plotting week 4 scatterplot...")
  res@xyMarkerColor     = "orange"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, extWaccm(week4,k), extOsiris(week4,k), res))       ; scatterplot
  
  print("Plotting week 3 scatterplot...")
  res@xyMarkerColor     = "yellow"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, extWaccm(week3,k), extOsiris(week3,k), res))       ; scatterplot
  
  print("Plotting week 2 scatterplot...")
  res@xyMarkerColor     = "green"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, extWaccm(week2,k), extOsiris(week2,k), res))       ; scatterplot
  
  print("Plotting week 1 scatterplot...")
  res@xyMarkerColor     = "blue"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, extWaccm(week1,k), extOsiris(week1,k), res))       ; scatterplot
  
  print("Plotting pre-Nabro scatterplot...")
  res@xyMarkerColor     = "gray"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, extWaccm(preNabro,k), extOsiris(preNabro,k), res))       ; scatterplot
  
  res@xyMarkLineModes   = "Lines"                ; choose which have markers
  res@xyMarkerColor     = "black"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, (/0,max(data)/), (/0,max(data)/), res))       ; scatterplot
  
  
  draw(plot)
  frame(wks)
  
  delete(wks)
  delete(plot)
  delete(res)
  
  ; Create a log scale for the bins.
  bin = fspan(0., 0.0005, 101)
;  bmin = 500.
;  bmax = 5e-6
;  nBins = 1000
;  dbin = (log(bmax) - log(bmin)) / (nBins - 1)
  
;  bin = fspan(0, nBins-1, nBins)
;  bin(:) = bmin*exp(bin(:) * dbin)
  
  pdffile=pdfpath+"hist_gsn_ext_"+sprintf("%04.1f",lev(k)/1000.)+"km"
  wks  = gsn_open_wks("pdf", pdffile)    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow")  ; choose a colormap

  res                                 = True
  res@tiMainString                    = "Extinction (750 nm) at "+sprintf("%4.1f",lev(k)/1000.)+" km"
  res@gsnHistogramSelectNiceIntervals = True   ; intervals now float
  res@gsnHistogramComputePercentages  = True    ; change left axis to %
  res@tmXBLabelAngleF                 = 270.    ; change label angle
  res@gsnHistogramCompare             = True
  res@gsnHistogramBinIntervals        = bin
  res@gsFillColor                     = "blue"
  res@tiXAxisString                   = "Extinction (/km)"

  print("Plotting gsn_histogram...")
  plot = gsn_histogram(wks, data, res)

  ; Now plot this histogram data in another way.
  numInBins = plot@NumInBins
  nBins = dimsizes(bin) - 1

  ; Normalize the data.
  totalCount = dim_sum(numInBins)
  pdf = int2flt(numInBins)
  pdf(0, :) = pdf(0, :) / totalCount(0) 
  pdf(1, :) = pdf(1, :) / totalCount(1) 

  ; Make the data step to how the bins.
  pdf2 = new((/ 2, nBins*2 /), "float")

  do i = 0, nBins - 1
    pdf2(:, 2*i)   = pdf(:, i)
    pdf2(:, 2*i+1) = pdf(:, i)
  end do

  xs = new(nBins*2, "float")
  xs(0) = bin(0)
  do i = 1, nBins-2
    xs(2*i-1) = bin(i)
    xs(2*i) = bin(i)
  end do
  xs(nBins*2-1) = bin(nBins-1)

  delete(wks)
  delete(plot)
  delete(res)  
  
  pdffile=pdfpath+"hist_ext_"+sprintf("%04.1f",lev(k)/1000.)+"km"
  wks  = gsn_open_wks("pdf", pdffile)    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow")  ; choose a colormap

  res                   = True
  res@gsnDraw           = False
  res@gsnFrame          = False
  res@tiMainString      = "Extinction (750 nm) at "+sprintf("%4.1f",lev(k)/1000.)+" km"
  res@xyLineThicknessF  = (/ 4.0, 4.0 /)
  res@xyLineColors      = (/ "black", "blue" /)
  res@xyDashPatterns    = (/ 0, 0 /)
  res@tiYAxisString     = "Probability (%)"
  res@tiXAxisString     = "Extinction (/km)"
;  res@trXLog            = True
;  res@trXMinF           = 5e-6
;  res@trXMaxF           = 500.
;  res@trYMaxF           = 15
;  res@tmXBMinorPerMajor = 8
  
  ; add a legend
  res@pmLegendDisplayMode    = "Always"              ; turn on legend
  res@pmLegendSide           = "Top"                 ; Change location of 
  res@pmLegendParallelPosF   = .225                   ; move units right
  res@pmLegendOrthogonalPosF = -0.3                  ; more neg = down
  res@pmLegendWidthF         = 0.12                  ; Change width and
  res@pmLegendHeightF        = 0.1                  ; height of legend.
  res@lgLabelFontHeightF     = .015                   ; change font height
  res@lgPerimOn              = False                 ; no box around
  res@xyExplicitLegendLabels = (/ " WACCM/CARMA", " Osiris" /)

  
  ave = new(2, "float")
  ave(0) = dim_avg(data(0,:))
  ave(1) = dim_avg(data(1,:))
  
  print(ave)
  
  print("Plotting histogram...")
  plot = gsn_csm_xy(wks, xs, pdf2 * 100, res)       ; contour the variable
  
  ; Draw vertical lines for the mean.
  res2                   = True
  res2@gsnDraw           = False
  res2@gsnFrame          = False
  res2@xyLineThicknessF  = (/ 2.5 /)
  res2@xyLineColors      = (/ "black" /)
  res2@xyDashPatterns    = (/ 2 /)
  
  overlay(plot, gsn_csm_xy(wks, (/ ave(0), ave(0) /), (/ 0, 100. /), res2))

  res2@xyLineColors      = (/ "blue" /)
  overlay(plot, gsn_csm_xy(wks, (/ ave(1), ave(1) /), (/ 0, 100. /), res2))

  draw(plot)
  frame(wks)

end

;correlateWACCMtoOsiris(16) ; 16.5 km
;correlateWACCMtoOsiris(17) ; 17.5 km
;correlateWACCMtoOsiris(18) ; 18.5 km
;correlateWACCMtoOsiris(19) ; 19.5 km
correlateWACCMtoOsiris(20) ; 20.5 km
correlateWACCMtoOsiris(21) ; 21.5 km
;correlateWACCMtoOsiris(22) ; 22.5 km
;correlateWACCMtoOsiris(23) ; 23.5 km
;correlateWACCMtoOsiris(24) ; 24.5 km
