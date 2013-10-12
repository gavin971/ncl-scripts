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
    
    if (exists .and. .not.(opt@overwritemass)) then
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
        dmass(:,:)=0.0e0
        dmass&time=time

        sulfur=new((/nTracers*nRuns,nTimes/),float)
        sulfur!1 = "time"
        sulfur@units = "kg"
        sulfur(:,:)=0.0e0
        sulfur&time=time
        printVarSummary(sulfur)        
      end if
      
      ntmin=min((/nt,nTimes/))-1
      t1=opt@timeoffset(i)
      t2=opt@timeoffset(i)+ntmin
      t3=opt@timeoffset2(i)
      print("t1="+t1+" t2="+t2+" t3="+t3+" ntmin="+ntmin)
      do j = 0, nTracers - 1
        k=j+i*nTracers
        tracer = mfileid->$tracerNames(j)$
        sulfur(k, t3:ntmin+t3) = tracer(t1:t2)*32.064999/tracerNames@molwts(j)
        if (j .ne. 0) then
          sulfur(k,:)=sulfur(k-1,:)+sulfur(k,:)
        end if
      end do
      so4_tot = mfileid->so4_tot
      dmass(i, t3:ntmin+t3) = so4_tot(t1:t2)
      delete(mfileid)
      delete(tracer)
    else
      print("mkdir -p "+massfilepath)
      a = systemfunc("mkdir -p "+massfilepath)
      historyPath = basePath + runNames(i) + "/atm/hist/*." + historyPattern + ".*"
      print("calling calculateMass, historyPath="+historyPath)
      test = calculateMass(historyPath, tracerNames, opt)
      print("return from calculateMass")
      printVarSummary(test)
      time=test&time
      nt=dimsizes(time)
      
      if (i .eq. 0) then
        nTimes=nt
        print("nTimes="+nTimes)

        dmass = new((/ nRuns, nTimes /), "float")
        dmass!0 = "run"
        dmass&run = runNames
        dmass!1 = "time"
        dmass@units = "kg"
        dmass(:,:)=0.0e0
        dmass&time=time

        sulfur=new((/nTracers*nRuns,nTimes/),float)
        sulfur!1 = "time"
        sulfur@units = "kg"
        sulfur(:,:)=0.0e0
        sulfur&time=time
        printVarSummary(sulfur)        
      end if
      
      print(" time@units:"+time@units+" nt="+nt)
      print(" dimsizes(test):"+dimsizes(test))
      ntmin=min((/nt,nTimes/))-1
      t1=opt@timeoffset(i)
      t2=opt@timeoffset(i)+ntmin
      t3=opt@timeoffset2(i)
      print("t1="+t1+" t2="+t2+" t3="+t3+" ntmin="+ntmin)

      print("rm "+massfilepath+massfile)
      rm=systemfunc("rm "+massfilepath+massfile)
      print("creating "+massfilepath+massfile)
      mfileid = addfile(massfilepath+massfile, "c")
      ; Create an UNLIMITED record dimension in the output netCDF file.  This is critical if 
      ;  the user plans to ever use NCO to concatenate the file along the time/record dimension.
;      filedimdef(mfileid,"time",-1,True)
      
      so4_tot = new(nt,float)
      so4_tot(:)=0.0e0
      do j = 0, nTracers - 1
        k=j+i*nTracers
        if (isStrSubset(tracerNames(j),"so4_a")) then
          so4_tot(:) = so4_tot(:) + (/test(j, :)/) ;- mass(i, j, 0)
        end if
        trcmass=test(j, :)
        trcmass!0="time"
        trcmass&time=time
        mfileid->$tracerNames(j)$=trcmass
      end do
      so4_tot!0="time"
      so4_tot&time=time
      mfileid->so4_tot=so4_tot
      delete(mfileid)

      printVarSummary(dmass)
      printVarSummary(so4_tot)
      print("i="+i+" t3="+t3+" ntmin="+ntmin+" t1="+t1+" t2="+t2)
      dmass(i,t3:ntmin+t3)=so4_tot(t1:t2)
      
      do j = 0, nTracers - 1
        k=j+i*nTracers
        print("k="+k+" j="+j+" "+tracerNames(j))
        sulfur(k,t3:ntmin+t3)=test(j,t1:t2)*32.064999/tracerNames@molwts(j)
        if (j .ne. 0) then
          sulfur(k,:)=sulfur(k-1,:)+sulfur(k,:)
        end if
      end do
      
      delete(test)
      delete(trcmass)
      
    end if
    
;    pdffile="dmass_"+runNames(i)+"_"+historyPattern
;    wks = gsn_open_wks("pdf", "$IMAGE_PATH/"+pdffile)
;    res = True
;    res@gsnMaximize           = True
;    res@vpHeightF             = 0.5              ; change aspect ratio
;    res@vpWidthF              = 1.0
;    res@xyMonoLineThickness = False
;    res@xyLineThicknesses = (/6.0,5.0,4.0,2.0,1.0/)
;    res@gsnLeftString = " "
;    res@gsnCenterString = "Mass burden (Tg)"
;    res@gsnRightString  = " "
;    res@xyLineColors    = (/"black","blue","red","green","orange"/)
;    res@xyMonoLineColor = False
;    res@xyMonoDashPattern = True
;    res@trXMaxF = 2.0
;    plot = gsn_csm_xy(wks, time/365., sulfur(:, :) /1e9, res)

;    delete(wks)
      
    delete(so4_tot)
    delete(time)
;    delete(res)
      
  end do ;i = 0, nRuns - 1
  
  dmass = dmass * 32.064999/115.108999
  ; Create a panel plot of the mass difference.
  print("")
  print(" Plotting ...")
  pdffile="dmass_Pinatubo_"+historyPattern
  wks = gsn_open_wks("pdf", "$IMAGE_PATH/"+pdffile)

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
  res@trXMinF = 1991.0
  res@trXMaxF = 1993.0
  res@tmXBMinorOn = False
;  res@tmXBMinorPerMajor = 11
  res@tmXBMaxTicks = 25
  res@tmXBMode= "Explicit"
  res@tmXBValues =fspan(1991,1993,25)
  res@tmXBLabelFontHeightF = 0.015
  res@tmXBLabels = (/"1991","F","M","A","M","J","J","A","S","O","N","D", \
                     "1992","F","M","A","M","J","J","A","S","O","N","D", \
                     "1993"/)
  plot = gsn_csm_xy(wks, 1991.+sulfur&time/365., sulfur(0:nTracers-1, :) /1e9, res)
  
  plot = gsn_csm_xy(wks, 1991.+sulfur&time/365., sulfur(:, :) /1e9, res)

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
  res@tmXBMaxTicks = 25
  res@tmXBMode= "Explicit"
  res@tmXBValues =fspan(1991,1993,25)
  res@tmXBLabelFontHeightF = 0.015
  res@tmXBLabels = (/"1991","F","M","A","M","J","J","A","S","O","N","D", \
                     "1992","F","M","A","M","J","J","A","S","O","N","D", \
                     "1993"/)
  plot = gsn_csm_xy(wks, 1991.+dmass&time/365., dmass(:, :) /1e9, res)

  delete(wks)

end


rootpath="/glade/scratch/mmills/archive/"
;rootpath="/glade/p/cesm/wawg/mmills/run/"
;caseid = (/ "b.e11.B55TRW5CN.f19_g16.Pin16-20.dgnh400nm.noqbo.001" /)
;caseid = (/ "b.e11.B55TRW5CN.f19_g16.Pin16-20.dgnh400nm.noqbo.001", \
;            "b.e11.B55TRW5CN.f19_g16.Pin16-20.dgnh1um.noqbo.001" , \
;            "b.e11.B55TRW5CN.f19_g16.Pin16-20.dgnh10um.noqbo.001" /)
;caseid="b.e11.B55TRW5CN.f19_g16.Pin16-18.dgnumhi1um.001"
;caseid = (/ "b.e11.B55TRW5CN.f19_g16.Pin16-18.dgnumhi1um.001", \
;            "b.e11.B55TRW5CN.f19_g16.Pin20-30.dgnh1um.noqbo.001" , \
;            "b.e11.B55TRW5CN.f19_g16.Pin17-25.dgnh1um.noqbo.001" /)
;         "b.e11.B55TRW5CN.f19_g16.Pin16-18.dgnh1um.qboon.001" /);, \
;         "b.e11.B55TRW5CN.f19_g16.MJMtest.013", \
;         "b.e11.B55TRW5CN.f19_g16.MJMtest.010" /)
;caseid = (/ "b.e12.B55TRW5CN.f19_g16.course.001", \
;            "b.e12.B55TRW5CN.f19_g16.CoursePin.002" /)
;caseid = (/ "b.e12.B55TRW5CN.f19_g16.CoursePin16-18.001", \
;            "b.e12.B55TRW5CN.f19_g16.CoursePin16-18.002" /)
caseid = (/ "b.e12.B55TRW5CN.f19_g16.CoursePin17-25.QBOon.001", \
            "b.e12.B55TRW5CN.f19_g16.CoursePin17-25.QBOon.002" /)
tracers=(/"SO2","H2SO4","so4_a2","so4_a1","so4_a3"/)
tracers@molwts=(/ 115.108999, 115.108999, 64.063799, 115.108999, 98.078479 /)
opt=False
opt@overwritemass=False
;opt@lev1=100.0
opt@lev1=1100.0
opt@timeoffset=(/0,0,0/)
opt@timeoffset2=(/0,0,0/)
;opt@timeoffset=(/0,0,0,3*12,12/)
;opt@timeoffset2=(/0,15,0,3*12,12/)
;plotMass( rootpath, "h0zm", caseid, tracers, opt )
plotMass( rootpath, "h3", caseid, tracers, opt )

;opt@overwritemass=True
;caseid="b.e11.B55TRW5CN.f19_g16.Pin20-30.dgnh1um.noqbo.001"
;plotMass( rootpath, "h0zm", caseid, tracers, opt )

;caseid="b.e11.B55TRW5CN.f19_g16.Pin16-18.dgnh1um.qboon.001"
;plotMass( rootpath, "h0zm", caseid, tracers, opt )

;opt@overwritemass=False
;caseid="b.e11.B55TRW5CN.f19_g16.MJMtest.013"
;plotMass( rootpath, "h0zm", caseid, tracers, opt )

;caseid="b.e11.B55TRW5CN.f19_g16.MJMtest.010"
;plotMass( rootpath, "h0zm", caseid, tracers, opt )
