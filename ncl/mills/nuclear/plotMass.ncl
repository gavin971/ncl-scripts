; Advection Test Case #1
; CAM-offline 3.1.3
;

load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "./calculateMass.ncl"


procedure plotMass(basePath:string, historyPattern:string, runNames[*]:string, \
  tracerNames[*]:string, opt:logical)
  
; opt: optional argument, true if IC file is to be added via opt@icfile
  
local nRuns, nTracers, historyPath, history, nTimes, mass, wks, res, resP, i

begin

  nRuns = dimsizes(runNames)
  nTracers = dimsizes(tracerNames)
  
  ; To figure out the size of the array that we will need, open the first
  ; set f history files.
  historyPath = basePath + runNames(0) + "/" + historyPattern + "/"
  print(historyPath)
  historyNames = systemfunc("ls " + historyPath)
;  print(historyNames)
  history = addfiles(historyPath+historyNames, "r")
  ListSetType(history, "cat")
;  nTimes = dimsizes(history[:]&time) 
  time1 = addfiles_GetVar(history, historyPath+historyNames, "time")
  nTimes = dimsizes(time1) 
  if (opt) then ; add time 0 for IC
    nTimes=nTimes+1
  end if
  print("nTimes="+nTimes)

  ; Create arrays to store all the mass data.
  mass = new((/ nRuns, nTracers, nTimes /), "float")
  mass!0 = "run"
  mass&run = runNames
  mass!1 = "tracer"
  mass&tracer = tracerNames
  mass!2 = "time"
  mass@units = "kg"
  
  dmass = new((/ nRuns, nTracers, nTimes /), "float")
  dmass!0 = "run"
  dmass&run = runNames
  dmass!1 = "tracer"
  dmass&tracer = tracerNames
  dmass!2 = "time"
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
    if (exists) then
      print("reading "+massfilepath+massfile)
      mfileid = addfile(massfilepath+massfile, "r")
      do j = 0, nTracers - 1
        dmass(i, j, :) = mfileid->$tracerNames(j)$
      end do
      delete(mfileid)
    else
      historyPath = basePath + runNames(i) + "/" + historyPattern + "/"
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

      print("creating "+massfilepath+massfile)
      mfileid = addfile(massfilepath+massfile, "c")
      do j = 0, nTracers - 1
        dmass(i, j, :) = (/mass(i, j, :)/) ;- mass(i, j, 0)

        print("  " + tracerNames(j) + " :  " + mass(i, j, nTimes-1) + " - " + \
          mass(i, j, 0) + " = " + dmass(i, j, nTimes-1))
        trcmass=dmass(i, j, :)
        trcmass!0="time"
        trcmass&time=time
        mfileid->$tracerNames(j)$=trcmass
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
  res@gsnFrame = False
  
  res@trYMinF = 0
  res@trYMaxF = 5.0
  
  res@trXMinF = 0
  res@trXMaxF = 10.0
  res@xyLineThicknessF = 6.0
  
;  res@tmYLFormat = "0*+^e"
;  res@tmYLPrecision = 3
  robock = (/ (/0.353, 4.948/), (/0.409, 4.892/), (/0.465, 4.837/), (/0.539, 4.795/), (/0.613, 4.768/), \
              (/0.687, 4.754/), (/0.761, 4.740/), (/0.873, 4.698/), (/0.984, 4.657/), (/1.188, 4.587/), \
              (/1.410, 4.504/), (/1.707, 4.421/), (/2.115, 4.310/), (/2.486, 4.129/), (/2.764, 4.005/), \
              (/3.060, 3.880/), (/3.468, 3.616/), (/3.783, 3.422/), (/4.099, 3.283/), (/4.358, 3.075/), \
              (/4.803, 2.798/), (/5.063, 2.673/), (/5.359, 2.493/), (/5.749, 2.243/), (/6.027, 2.118/), \
              (/6.305, 1.938/), (/6.583, 1.771/), (/6.787, 1.660/), (/7.139, 1.535/), (/7.510, 1.341/), \
              (/7.862, 1.202/), (/8.103, 1.119/), (/8.363, 0.994/), (/8.789, 0.856/), (/9.067, 0.772/), \
              (/9.364, 0.689/), (/9.568, 0.634/), (/9.791, 0.578/), (/9.939, 0.550/), (/10.050, 0.537/) /)
  
  robock(:,0)=robock(:,0)-134./365. ; Robock starts May 15, vs Jan 1

  
  do i = 0, nRuns - 1
    res@gsnLeftString = " "
    res@gsnCenterString = "Mass burden (Tg)"
    res@gsnRightString = " "
;    plot(i) = gsn_csm_xy(wks, dmass&time, dmass(i, nTracers - 1, :) / mass(i, nTracers-1, 0), res)
    res@xyLineColor       = "black"
    plot = gsn_csm_xy(wks, dmass&time/365., dmass(i, nTracers - 1, :) /1e9, res)
    res@xyLineColor       = "red"
;    print(robock(:,0))
;    print(robock(:,1))
    plot = gsn_csm_xy(wks, robock(:,0), robock(:,1), res)
    frame(wks)
  end do
  
;  resP = True
;  resP@txString = "Total Mass"
;  gsn_panel(wks, plot, (/ nRuns, 1 /), resP)
  delete(wks)
  
  host   = systemfunc("hostname")
  ostype = systemfunc("echo $OSTYPE")
  if (host.eq."modeling1") then
    muttcmd="mutt mikejmills@mac.com -a $IMAGE_PATH/"+pdffile+".pdf -s "+pdffile+" < /dev/null"
    print(muttcmd)
    mutter=systemfunc(muttcmd)
    print(mutter)
  end if
  
end

rootpath="$MODEL_PATH/cesm/case/"
opt=True
opt@icfile=rootpath+"/b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg/ic/b40.rcp4_5.2deg.wcm.carma.bc5tg.IndPak.001.cam2.i.2013-01-01-00000.nc"
plotMass( rootpath, "h0zm", "b40.rcp4_5.2deg.wcm.carma.bc5tgIP.003", "SOOT01", opt )
