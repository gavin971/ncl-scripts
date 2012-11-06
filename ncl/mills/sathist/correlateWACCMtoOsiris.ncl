; Comparison of WACCM Nabro extinction calculations with Osiris observations
; Set the following environment variables
; hspath : full path to the input sathist file
; outfile : file created by interWACCMtoOsiris.ncl, i.e. $case.cam.hs.OsirisLevs.nc
; caseid : used for the name of the PDF file generated

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("correlateWACCMtoOsiris")
procedure correlateWACCMtoOsiris()

begin

  print("")
  k1=getenv("k")
  
  if ismissing(k1) then
    k = 19
    print("default k="+k+", reading data...")
  else 
;    print(typeof(k))
    k=stringtointeger(k1)
    print("k="+k+", reading data...")
  end if

  ; Get the osiris and WACCM data
;  waccmRun="f40.Nabro.2deg.wcm.sd.carma.sulf.006"
  waccmRun=getenv("caseid")
  print("waccmRun="+waccmRun)
  
  waccmRun=getenv("caseid")
  print("waccmRun="+waccmRun)
  
;  pathList=(/"/Volumes/MillsPassport/cesm/case/",\
;             "/Volumes/Data/Models/cesm/case/"/)

;  i=0
;  filePresent=False           
;  do while(i.lt.dimsizes(pathList).and..not.(filePresent))
;    rootdir=pathList(i)+waccmRun+"/"
;    filePresent=isfilepresent(rootdir)
;    print(i+" "+rootdir+" "+filePresent)
;    i=i+1
;  end do

;  if (.not.(filePresent)) then
;    print("file not found")
;    exit
;  end if
  
  hspath=getenv("runpath")
  hspath=hspath+"/"+waccmRun+"/hs"
  print("hspath="+hspath)
  waccmFile=waccmRun+".cam.hs.OsirisLevs.nc"
  print("waccmFile="+waccmFile)
  
  osirisfile="$DATA_PATH/Osiris/osiris_hs_nabro_20110601_20110822.nc"
  print("adding "+osirisfile)
  osiris=addfile(osirisfile,"r")
;  osirisFile = addfile("$DATA_PATH/Osiris/osiris_hs_nabro_20110601_20110822.nc","r")
  extOsiris  = osiris->AEXT
  exterrOsiris  = osiris->AEXTerr
  
;  waccmFile = addfile(rootdir+"hs/"+waccmRun+".cam.hs.OsirisLevs.nc","r")
  hsfile=hspath+"/"+waccmFile
  print("adding "+hsfile)
  waccm = addfile(hsfile,"r")
  extWaccm  = waccm->PURSULEX750

  dimsOsiris=dimsizes(extOsiris)
  dims=dimsizes(extWaccm)
  ncol=min((/dims(0),dimsOsiris(0)/))
  nlev=dims(1)
  print("extOsiris:"+dimsOsiris)
  print("extWaccm:"+dims)
  print(" ")
  
  lev = osiris->lev
  date = osiris->date
  date_waccm = waccm->date
  prof_num_osiris=osiris->prof_num
  
  prof_num_waccm=waccm->prof_num
  nprofsWaccm=dimsizes(prof_num_waccm)
  print("prof_num_waccm: "+nprofsWaccm+" "+prof_num_waccm(0)+" "+prof_num_waccm(nprofsWaccm-1))
  print("date_waccm: "+dimsizes(date_waccm)+" "+date_waccm(0)+" "+date_waccm(nprofsWaccm-1))
  print(" ")
  
  preNabroOsiris=ind(prof_num_osiris.ge.prof_num_waccm(0).and.date.lt.20110613)
  nprofsPreNabroOsiris=dimsizes(preNabroOsiris)
  print("preNabro Osiris:"+nprofsPreNabroOsiris+" "+prof_num_osiris(preNabroOsiris(0))+" "+prof_num_osiris(preNabroOsiris(0)+nprofsPreNabroOsiris-1))

  preNabroWaccm=ind(prof_num_waccm.le.prof_num_osiris(preNabroOsiris(0)+nprofsPreNabroOsiris-1))
  nprofsPreNabroWaccm=dimsizes(preNabroWaccm)
  print("preNabro Waccm:"+nprofsPreNabroWaccm+" "+prof_num_waccm(preNabroWaccm(0))+" "+prof_num_waccm(preNabroWaccm(0)+nprofsPreNabroWaccm-1))
  print(" ")

  if ((nprofsPreNabroOsiris.ne.nprofsPreNabroWaccm).or.\
      (prof_num_osiris(preNabroOsiris(0)).ne.prof_num_waccm(preNabroWaccm(0))).or.\
      (prof_num_osiris(preNabroOsiris(0)+nprofsPreNabroOsiris-1).ne.prof_num_waccm(preNabroWaccm(0)+nprofsPreNabroWaccm-1))) then      
    print("Misalignment! Aborting.")
    return    
  end if

  Week1Osiris=ind(date.ge.20110613.and.date.lt.20110620)
  nprofsWeek1Osiris=dimsizes(Week1Osiris)
  print("Week 1 Osiris:"+nprofsWeek1Osiris+" "+prof_num_osiris(Week1Osiris(0))+" "+prof_num_osiris(Week1Osiris(0)+nprofsWeek1Osiris-1))
  Week1Waccm=ind(prof_num_waccm.ge.prof_num_osiris(Week1Osiris(0)).and.prof_num_waccm.le.prof_num_osiris(Week1Osiris(0)+nprofsWeek1Osiris-1))
  nprofsWeek1Waccm=dimsizes(Week1Waccm)
  print("Week 1 Waccm:"+nprofsWeek1Waccm+" "+prof_num_waccm(Week1Waccm(0))+" "+prof_num_waccm(Week1Waccm(0)+nprofsWeek1Waccm-1))
  print(" ")

  if ((nprofsWeek1Osiris.ne.nprofsWeek1Waccm).or.\
      (prof_num_osiris(Week1Osiris(0)).ne.prof_num_waccm(Week1Waccm(0))).or.\
      (prof_num_osiris(Week1Osiris(0)+nprofsWeek1Osiris-1).ne.prof_num_waccm(Week1Waccm(0)+nprofsWeek1Waccm-1))) then      
    print("Misalignment! Aborting.")
    return    
  end if
  
  Week2Osiris=ind(date.ge.20110620.and.date.lt.20110627)
  nprofsWeek2Osiris=dimsizes(Week2Osiris)
  print("Week 2 Osiris:"+nprofsWeek2Osiris+" "+prof_num_osiris(Week2Osiris(0))+" "+prof_num_osiris(Week2Osiris(0)+nprofsWeek2Osiris-1))
  Week2Waccm=ind(prof_num_waccm.ge.prof_num_osiris(Week2Osiris(0)).and.prof_num_waccm.le.prof_num_osiris(Week2Osiris(0)+nprofsWeek2Osiris-1))
  nprofsWeek2Waccm=dimsizes(Week2Waccm)
  print("Week 2 Waccm:"+nprofsWeek2Waccm+" "+prof_num_waccm(Week2Waccm(0))+" "+prof_num_waccm(Week2Waccm(0)+nprofsWeek2Waccm-1))
  print(" ")

  if ((nprofsWeek2Osiris.ne.nprofsWeek2Waccm).or.\
      (prof_num_osiris(Week2Osiris(0)).ne.prof_num_waccm(Week2Waccm(0))).or.\
      (prof_num_osiris(Week2Osiris(0)+nprofsWeek2Osiris-1).ne.prof_num_waccm(Week2Waccm(0)+nprofsWeek2Waccm-1))) then      
    print("Misalignment! Aborting.")
    return    
  end if
  
  Week3Osiris=ind(date.ge.20110627.and.date.lt.20110704)
  nprofsWeek3Osiris=dimsizes(Week3Osiris)
  print("week 3 Osiris:"+nprofsWeek3Osiris+" "+prof_num_osiris(Week3Osiris(0))+" "+prof_num_osiris(Week3Osiris(0)+nprofsWeek3Osiris-1))
  Week3Waccm=ind(prof_num_waccm.ge.prof_num_osiris(Week3Osiris(0)).and.prof_num_waccm.le.prof_num_osiris(Week3Osiris(0)+nprofsWeek3Osiris-1))
  nprofsWeek3Waccm=dimsizes(Week3Waccm)
  print("week 3 Waccm:"+nprofsWeek3Waccm+" "+prof_num_waccm(Week3Waccm(0))+" "+prof_num_waccm(Week3Waccm(0)+nprofsWeek3Waccm-1))
  print(" ")

  if ((nprofsWeek3Osiris.ne.nprofsWeek3Waccm).or.\
      (prof_num_osiris(Week3Osiris(0)).ne.prof_num_waccm(Week3Waccm(0))).or.\
      (prof_num_osiris(Week3Osiris(0)+nprofsWeek3Osiris-1).ne.prof_num_waccm(Week3Waccm(0)+nprofsWeek3Waccm-1))) then      
    print("Misalignment! Aborting.")
    return    
  end if
  
  Week4Osiris=ind(date.ge.20110704.and.date.lt.20110711)
  nprofsWeek4Osiris=dimsizes(Week4Osiris)
  print("week 4 Osiris:"+nprofsWeek4Osiris+" "+prof_num_osiris(Week4Osiris(0))+" "+prof_num_osiris(Week4Osiris(0)+nprofsWeek4Osiris-1))
  Week4Waccm=ind(prof_num_waccm.ge.prof_num_osiris(Week4Osiris(0)).and.prof_num_waccm.le.prof_num_osiris(Week4Osiris(0)+nprofsWeek4Osiris-1))
  nprofsWeek4Waccm=dimsizes(Week4Waccm)
  print("week 4 Waccm:"+nprofsWeek4Waccm+" "+prof_num_waccm(Week4Waccm(0))+" "+prof_num_waccm(Week4Waccm(0)+nprofsWeek4Waccm-1))
  print(" ")

  if ((nprofsWeek4Osiris.ne.nprofsWeek4Waccm).or.\
      (prof_num_osiris(Week4Osiris(0)).ne.prof_num_waccm(Week4Waccm(0))).or.\
      (prof_num_osiris(Week4Osiris(0)+nprofsWeek4Osiris-1).ne.prof_num_waccm(Week4Waccm(0)+nprofsWeek4Waccm-1))) then      
    print("Misalignment! Aborting.")
    return    
  end if
  
  PostWeek4Osiris=ind(date.ge.20110711)
  nprofsPostWeek4Osiris=dimsizes(PostWeek4Osiris)
  print("post week 4 Osiris:"+nprofsPostWeek4Osiris+" "+prof_num_osiris(PostWeek4Osiris(0))+" "+prof_num_osiris(PostWeek4Osiris(0)+nprofsPostWeek4Osiris-1))
  PostWeek4Waccm=ind(prof_num_waccm.ge.prof_num_osiris(PostWeek4Osiris(0)).and.prof_num_waccm.le.prof_num_osiris(PostWeek4Osiris(0)+nprofsPostWeek4Osiris-1))
  nprofsPostWeek4Waccm=dimsizes(PostWeek4Waccm)
  print("post week 4 Waccm:"+nprofsPostWeek4Waccm+" "+prof_num_waccm(PostWeek4Waccm(0))+" "+prof_num_waccm(PostWeek4Waccm(0)+nprofsPostWeek4Waccm-1))
  print(" ")

  if ((nprofsPostWeek4Osiris.ne.nprofsPostWeek4Waccm).or.\
      (prof_num_osiris(PostWeek4Osiris(0)).ne.prof_num_waccm(PostWeek4Waccm(0))).or.\
      (prof_num_osiris(PostWeek4Osiris(0)+nprofsPostWeek4Osiris-1).ne.prof_num_waccm(PostWeek4Waccm(0)+nprofsPostWeek4Waccm-1))) then      
    print("Misalignment! Aborting.")
    return    
  end if
    
  extOsiris  = where(ismissing(extOsiris),0.,extOsiris)
    
  data = new((/ 2, ncol /), "float")
  data(0, :) = extWaccm(0:ncol-1,k)
  data(1, :) = (/extOsiris(0:ncol-1,k)/)
  data@_FillValue = -999.
  
  ; Plot the correlation of the data.
  pdfpath="$IMG_PATH/"+waccmRun+"/Osiris/"
  a = systemfunc("mkdir -p "+pdfpath)
  pdffile=pdfpath+waccmRun+".corr_ext_"+sprintf("%04.1f",lev(k)/1000.)+"km"
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
  plot = gsn_csm_xy(wks, extWaccm(PostWeek4Waccm,k), extOsiris(PostWeek4Osiris,k), res)       ; scatterplot
  
  print("Plotting week 4 scatterplot...")
  res@xyMarkerColor     = "orange"                    ; Marker color
;  plot = gsn_csm_xy(wks, extWaccm(week4,k), extOsiris(week4,k), res)       ; scatterplot
  overlay(plot, gsn_csm_xy(wks, extWaccm(Week4Waccm,k), extOsiris(Week4Osiris,k), res))       ; scatterplot
  
  print("Plotting week 3 scatterplot...")
  res@xyMarkerColor     = "yellow"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, extWaccm(Week3Waccm,k), extOsiris(Week3Osiris,k), res))       ; scatterplot
  
  print("Plotting week 2 scatterplot...")
  res@xyMarkerColor     = "green"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, extWaccm(Week2Waccm,k), extOsiris(Week2Osiris,k), res))       ; scatterplot
  
  print("Plotting week 1 scatterplot...")
  res@xyMarkerColor     = "blue"                    ; Marker color
  res@xyMarkerSizeF     = 0.005                     ; Marker size (default 0.01)
  overlay(plot, gsn_csm_xy(wks, extWaccm(Week1Waccm,k), extOsiris(Week1Osiris,k), res))       ; scatterplot
  
  print("Plotting pre-Nabro scatterplot...")
  res@xyMarkerColor     = "gray"                    ; Marker color
  res@xyMarkerSizeF     = 0.002                     ; Marker size (default 0.01)
  overlay(plot, gsn_csm_xy(wks, extWaccm(preNabroWaccm,k), extOsiris(preNabroOsiris,k), res))       ; scatterplot
  
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
  
  pdffile=pdfpath+waccmRun+".hist_gsn_ext_"+sprintf("%04.1f",lev(k)/1000.)+"km"
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
  
  pdffile=pdfpath+waccmRun+".hist_ext_"+sprintf("%04.1f",lev(k)/1000.)+"km"
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
;correlateWACCMtoOsiris(20) ; 20.5 km
;correlateWACCMtoOsiris(21) ; 21.5 km
;correlateWACCMtoOsiris(22) ; 22.5 km
;correlateWACCMtoOsiris(23) ; 23.5 km
;correlateWACCMtoOsiris(24) ; 24.5 km

correlateWACCMtoOsiris()
