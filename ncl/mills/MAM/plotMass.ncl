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
  
  ; To figure out the size of the array that we will need, open the first
  ; set f history files.
  historyPath = basePath + runNames(0) + "/atm/hist/*." + historyPattern + ".*"
  print(historyPath)
  historyNames = systemfunc("ls " + historyPath)
;  print(historyNames)
;  history = addfiles(historyPath+historyNames, "r")
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
;  nTimes = dimsizes(history[:]&time) 
;  time1 = addfiles_GetVar(history, historyPath+historyNames, "time")
  time1 = addfiles_GetVar(history, historyNames, "time")
  nTimes = dimsizes(time1) 
  print("nTimes="+nTimes)
  
  dmass = new((/ nRuns, nTimes /), "float")
  dmass!0 = "run"
  dmass&run = runNames
  dmass!1 = "time"
  dmass@units = "kg"
  dmass(:,:)=0.0e0
  dmass&time=time1
  
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
      sulfur=new((/nTracers,nt/),float)
      do j = 0, nTracers - 1
        sulfur(j, :) = mfileid->$tracerNames(j)$
        sulfur(j, :) = sulfur(j, :)*32.064999/tracerNames@molwts(j)
        if (j .ne. 0) then
          sulfur(j,:)=sulfur(j-1,:)+sulfur(j,:)
        end if
      end do
      so4_tot = mfileid->so4_tot
      ntmin=min((/nt,nTimes/))-1
      dmass(i, 0:ntmin) = so4_tot(opt@timeoffset(i):opt@timeoffset(i)+ntmin)
      delete(mfileid)
    else
      print("mkdir -p "+massfilepath)
      a = systemfunc("mkdir -p "+massfilepath)
      historyPath = basePath + runNames(i) + "/atm/hist/*." + historyPattern + ".*"
      print("calling calculateMass, historyPath="+historyPath)
      test = calculateMass(historyPath, tracerNames, opt)
      sulfur=test
      sulfur(:,:)=0.0
      time=test&time
      nt=dimsizes(time)
      print(" time@units:"+time@units+" nt="+nt)
      print(" dimsizes(test):"+dimsizes(test))

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
        if (isStrSubset(tracerNames(j),"so4_a")) then
          so4_tot(:) = so4_tot(:) + (/test(j, :)/) ;- mass(i, j, 0)
        end if
        trcmass=test(j, :)
        trcmass!0="time"
        trcmass&time=time
        mfileid->$tracerNames(j)$=trcmass
        sulfur(j,:)=test(j,:)*32.064999/tracerNames@molwts(j)
        if (j .ne. 0) then
          sulfur(j,:)=sulfur(j-1,:)+sulfur(j,:)
        end if
      end do
      so4_tot!0="time"
      so4_tot&time=time
      mfileid->so4_tot=so4_tot
      ntmin=min((/nt,nTimes/))-1
      dmass(i,0:ntmin)=so4_tot(opt@timeoffset(i):opt@timeoffset(i)+ntmin)
      
      delete(mfileid)
      delete(test)
      delete(trcmass)
      
    end if
    
    pdffile="dmass_"+runNames(i)+"_"+historyPattern
    wks = gsn_open_wks("pdf", "$IMAGE_PATH/"+pdffile)
    res = True
    res@gsnMaximize           = True
    res@vpHeightF             = 0.5              ; change aspect ratio
    res@vpWidthF              = 1.0
    res@xyMonoLineThickness = False
    res@xyLineThicknesses = (/6.0,5.0,4.0,2.0,1.0/)
    res@gsnLeftString = " "
    res@gsnCenterString = "Mass burden (Tg)"
    res@gsnRightString  = " "
    res@xyLineColors    = (/"black","blue","red","green","orange"/)
    res@xyMonoLineColor = False
    res@xyMonoDashPattern = True
    res@trXMaxF = 2.0
    plot = gsn_csm_xy(wks, time/365., sulfur(:, :) /1e9, res)

    delete(wks)
      
    delete(so4_tot)
    delete(time)
    delete(sulfur)
    delete(res)
      
  end do
  
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

;  res@trYMinF = 0.0
;  res@trYMaxF = 13.0

;  res@trXMinF = 0
  res@xyLineColors    = (/"black","blue","red"/)
  res@xyMonoLineColor = False
  res@xyMonoDashPattern = True
  res@trXMaxF = 2.0
  res@xyLineThicknessF = 6.0

  res@gsnLeftString = " "
  res@gsnCenterString = "Mass burden (Tg)"
  res@gsnRightString = " "
  plot = gsn_csm_xy(wks, dmass&time/365., dmass(:, :) /1e9, res)

  delete(wks)

end


rootpath="/glade/scratch/mmills/archive/"
;caseid="b.e11.B55TRW5CN.f19_g16.Pin16-18.dgnumhi1um.001"
caseid=(/"b.e11.B55TRW5CN.f19_g16.Pin16-18.dgnumhi1um.001", \
         "b.e11.B55TRW5CN.f19_g16.Pin20-30.dgnh1um.noqbo.001", \
         "b.e11.B55TRW5CN.f19_g16.Pin16-18.dgnh1um.qboon.001" /);, \
;         "b.e11.B55TRW5CN.f19_g16.MJMtest.013", \
;         "b.e11.B55TRW5CN.f19_g16.MJMtest.010" /)
tracers=(/"H2SO4","SO2","so4_a1","so4_a2","so4_a3"/)
tracers@molwts=(/ 115.108999, 115.108999, 115.108999, 64.063799, 98.078479 /)
opt=False
opt@overwritemass=False
;opt@lev1=100.0
opt@lev1=1100.0
opt@timeoffset=(/0,0,0,3*12,12/)
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
