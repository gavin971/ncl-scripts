load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "../nuclear/calculateMass.ncl"

procedure plotMass(basePath:string, historyPattern:string, runNames[*]:string, \
  tracerNames[*]:string, opt:logical)
  
; opt: optional argument, true if IC file is to be added via opt@icfile
  
local nRuns, nTracers, historyPath, history, nTimes, mass, wks, res, resP, i

begin

  nRuns = dimsizes(runNames)
  nTracers = dimsizes(tracerNames)
    
  ; Calculate the tracer masses for all the runs.
  do i = 0, nRuns - 1
    print("")
    print("")
    print("Run :  " + runNames(i))
    
    massfilepath = basePath + runNames(i) + "/atm/proc/"
    massfile = runNames(i) + "." + historyPattern + ".mass.nc"    
    exists = isfilepresent(massfilepath+massfile)
    
    if (exists) then
      print("reading "+massfilepath+massfile)
      mfileid = addfile(massfilepath+massfile, "r")
      time=mfileid->time
      nt=dimsizes(time)
      print("  nt="+nt+" units="+time@units)
            
      if (i .eq. 0) then
        nTimes=nt
        print("nTimes="+nTimes)

        dmass = new((/ nRuns, nTimes /), "float")
        dmass!0 = "run"
        dmass&run = runNames
        dmass!1 = "time"
        dmass@units = "kg"
;        dmass(:,:)=0.0e0
        dmass&time=ispan(1,nTimes,1)

        sulfur=new((/nTracers*nRuns,nTimes/),float)
        sulfur!1 = "time"
        sulfur@units = "kg"
;        sulfur(:,:)=0.0e0
        sulfur&time=dmass&time
        printVarSummary(sulfur)        
      end if
      
      ntmin=nTimes-1
      t1=opt@timeoffset(i)

      print("t1="+t1+" ntmin="+ntmin)
      do j = 0, nTracers - 1
        k=j+i*nTracers
        tracer = mfileid->$tracerNames(j)$
        sulfur(k, t1:ntmin) = tracer(:)*32.064999/tracerNames@molwts(j)
        if (j .ne. 0) then
          sulfur(k,:)=sulfur(k-1,:)+sulfur(k,:)
        end if
      end do
      so4_tot = mfileid->so4_tot
      dmass(i, t1:ntmin) = so4_tot(:)
      delete(mfileid)
      delete(tracer)
    else
      print("file not found: "+massfilepath)
      exit
    end if
    
    delete(so4_tot)
    delete(time)
      
  end do ;i = 0, nRuns - 1
  
  dmass = dmass * 32.064999/115.108999
  ; Create a panel plot of the mass difference.
  print("")
  print(" Plotting ...")
  pdffile="dmass_GeoMIP_"+historyPattern
  pdfpath="$IMAGE_PATH/"+runNames(0)+"/"
  mkdir=systemfunc("mkdir -p "+pdfpath)
  wks = gsn_open_wks("pdf", pdfpath+pdffile)

  res = True
  res@gsnMaximize           = True
  res@vpHeightF             = 0.5              ; change aspect ratio
  res@vpWidthF              = 1.0
  res@xyMonoLineThickness = True
  res@xyLineThicknessF = 4.0
;  res@xyLineThicknesses = (/6.0,5.0,4.0,2.0,1.0,\
;                            6.0,5.0,4.0,2.0,1.0,\
;                            6.0,5.0,4.0,2.0,1.0/)
  res@gsnLeftString = " "
  res@gsnCenterString = " "
  res@gsnRightString  = " "
  res@tiYAxisString = "Mass burden (Tg S)"
  res@tiXAxisString = "Year"
  res@xyLineColors    = (/"orange","red","green4","blue","black",\
                          "orange","red","green4","blue","black",\
                          "orange","red","green4","blue","black"/)
  res@xyMonoLineColor = False
  res@xyMonoDashPattern = False
  res@xyDashPatterns = (/0,0,0,0,0,1,1,1,1,1,2,2,2,2,2/)
;  res@trXMinF = 1991.0
;  res@trXMaxF = 1993.0
  res@tmXBMinorOn = False
;  res@tmXBMinorPerMajor = 11
;  res@tmXBMaxTicks = 25
;  res@tmXBMode= "Explicit"
;  res@tmXBValues =fspan(1991,1993,25)
;  res@tmXBLabelFontHeightF = 0.015
;  tickLabels = (/"1991","F","M","A","M","J","J","A","S","O","N","D", \
;                     "1992","F","M","A","M","J","J","A","S","O","N","D", \
;                     "1993"/)
;  res@tmXBLabels = tickLabels
  plot = gsn_csm_xy(wks, 2012.+sulfur&time/365., sulfur(0:nTracers-1, :) /1e9, res)
  
  plot = gsn_csm_xy(wks, 2012.+sulfur&time/365., sulfur(:, :) /1e9, res)

  delete(res)
  res = True
  res@gsnMaximize           = True
  res@vpHeightF             = 0.5              ; change aspect ratio
  res@vpWidthF              = 1.0

;  res@trYMinF = 0.0
;  res@trYMaxF = 13.0

;  res@trXMinF = 0
;  res@xyLineColors    = (/"black","blue","red"/)
  res@xyLineColor    = "black"
  res@xyMonoLineColor = True
  res@xyMonoDashPattern = False
  res@xyDashPatterns = (/0,1,2/)
  res@trXMinF = 1991.0
  res@trXMaxF = 1993.0
  res@xyLineThicknessF = 4.0

  res@gsnLeftString = " "
  res@gsnCenterString = " "
  res@gsnRightString = " "
  res@tiYAxisString = "Sulfate mass burden (Tg S)"
  res@tiXAxisString = "Year"
  res@tmXBMinorOn = False
;  res@tmXBMinorPerMajor = 11
;  res@tmXBMaxTicks = 25
;  res@tmXBMode= "Explicit"
;  res@tmXBValues =fspan(1991,1993,25)
  res@tmXBLabelFontHeightF = 0.015
;  res@tmXBLabels = tickLabels
  plot = gsn_csm_xy(wks, 2012.+dmass&time/365., dmass(:, :) /1e9, res)

  delete(wks)

end


rootpath="/glade/scratch/mmills/archive/"
caseid = (/ \
       "b.e12.BRCP45W5CN.f19_g16.control.002", \
       "b.e12.BRCP45W5CN.f19_g16.geo01Tgpy2015TR.002", \
       "b.e12.BRCP45W5CN.f19_g16.geo02Tgpy2015TR.002", \
       "b.e12.BRCP45W5CN.f19_g16.geo04Tgpy2015TR.002", \
       "b.e12.BRCP45W5CN.f19_g16.geo10Tgpy2015TR.002" \
       /)
                     
tracers=(/"SO2","H2SO4","so4_a2","so4_a1","so4_a3"/)
tracers@molwts=(/ 64.063799, 98.078479, 115.108999, 115.108999, 115.108999 /)

opt=False
opt@timeoffset=(/0,36,36,36,36/)

plotMass( rootpath, "h0", caseid, tracers, opt )

