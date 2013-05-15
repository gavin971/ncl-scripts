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
  historyPath = basePath + runNames(0) + "/atm/proc/" + historyPattern + "/"
  print(historyPath)
  historyNames = systemfunc("ls " + historyPath)
;  print(historyNames)
  history = addfiles(historyPath+historyNames, "r")
  ListSetType(history, "cat")
;  nTimes = dimsizes(history[:]&time) 
  time1 = addfiles_GetVar(history, historyPath+historyNames, "time")
  nTimes = dimsizes(time1) 
  print("nTimes="+nTimes)

  ; Create arrays to store all the mass data.
  mass = new((/ nRuns, nTracers, nTimes /), "float")
  mass!0 = "run"
  mass&run = runNames
  mass!1 = "tracer"
  mass&tracer = tracerNames
  mass!2 = "time"
  mass@units = "kg"
  
  dmass = new((/ nRuns, nTimes /), "float")
  dmass!0 = "run"
  dmass&run = runNames
  dmass!1 = "time"
  dmass@units = "kg"
  
  ; Calculate the tracer masses for all the runs.
  do i = 0, nRuns - 1
    print("")
    print("")
    print("Run :  " + runNames(i))
    
    massfilepath = basePath + runNames(i) + "/proc/"
    massfile = runNames(i) +".mass.nc"
    print("mkdir -p "+massfilepath)
    a = systemfunc("mkdir -p "+massfilepath)
    
    exists = isfilepresent(massfilepath+massfile)
    if (exists .and. .not.(opt@overwritemass)) then
      print("reading "+massfilepath+massfile)
      mfileid = addfile(massfilepath+massfile, "r")
      do j = 0, nTracers - 1
        mass(i, j, :) = mfileid->$tracerNames(j)$
      end do
      dmass(i, :) = mfileid->so4_tot
      delete(mfileid)
    else
      historyPath = basePath + runNames(i) + "/atm/proc/" + historyPattern + "/"
      print("calling calculateMass, historyPath="+historyPath)
      test = calculateMass(historyPath, tracerNames, opt)
      time=test&time
      time@units="days"
      dmass&time = time
      mass&time = time
      print(dimsizes(test))
      print(dimsizes(mass))
      print(i)
      mass(i, :, :) = (/test/)

      print("rm "+massfilepath+massfile)
      rm=systemfunc("rm "+massfilepath+massfile)
      print("creating "+massfilepath+massfile)
      mfileid = addfile(massfilepath+massfile, "c")
      ; Create an UNLIMITED record dimension in the output netCDF file.  This is critical if 
      ;  the user plans to ever use NCO to concatenate the file along the time/record dimension.
;      filedimdef(mfileid,"time",-1,True)
      
      dmass(:,:)=0.0e0
      do j = 0, nTracers - 1
        dmass(i, :) = dmass(i, :) + (/mass(i, j, :)/) ;- mass(i, j, 0)

        print("  " + tracerNames(j) + " :  " + mass(i, j, nTimes-1) + " - " + \
          mass(i, j, 0) + " = " + mass(i, j, nTimes-1))
        trcmass=mass(i, j, :)
        trcmass!0="time"
        trcmass&time=time
        so4_tot=dmass(i,:)
        so4_tot!0="time"
        so4_tot&time=time
        mfileid->$tracerNames(j)$=trcmass
        mfileid->so4_tot=so4_tot
      end do
      
      delete(mfileid)
      
    end if
  end do
  
  ; Create a panel plot of the mass difference.
  print("")
  print(" Plotting ...")
  pdffile="dmass_"+runNames(0)+"_"+historyPattern
  wks = gsn_open_wks("pdf", "$IMAGE_PATH/"+pdffile)
;  plot = new(nRuns, graphic)

  res = True
  res@gsnMaximize           = True
  res@vpHeightF             = 0.5              ; change aspect ratio
  res@vpWidthF              = 1.0
;  res@gsnDraw = False
;  res@gsnFrame = False

  res@trYMinF = 0.19
  res@trYMaxF = 0.26

;  res@trXMinF = 0
;  res@trXMaxF = 10.0
  res@xyLineThicknessF = 6.0

  res@gsnLeftString = " "
  res@gsnCenterString = "Mass burden (Tg)"
  res@gsnRightString = " "
  res@xyLineColor       = "black"
  plot = gsn_csm_xy(wks, dmass&time/365., dmass(0, :) /1e9, res)

  delete(wks)

end


rootpath="/glade/scratch/mmills/archive/"
caseid="b.e11.B55TRW5CN.f19_g16.MJMtest.011"
tracers=(/"so4_a1","so4_a2","so4_a3"/)
opt=False
opt@overwritemass=False
opt@lev1=100.0
plotMass( rootpath, "h0zm", caseid, tracers, opt )

caseid="b.e11.B55TRW5CN.f19_g16.MJMtest.013b"
plotMass( rootpath, "h0zm", caseid, tracers, opt )
