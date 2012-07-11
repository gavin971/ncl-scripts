; Advection Test Case #1
; CAM-offline 3.1.3
;

load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "./calculateMass.ncl"


procedure plotMass(basePath:string, historyPattern:string, runNames[*]:string, \
  tracerNames[*]:string)
  
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
  nTimes = dimsizes(history[:]&time) 

  ; Create arrays to store all the mass data.
  mass = new((/ nRuns, nTracers, nTimes /), "float")
  mass!0 = "run"
  mass&run = runNames
  mass!1 = "tracer"
  mass&tracer = tracerNames
  mass!2 = "time"
  mass&time = history[:]&time
  mass@units = "kg"
  
  dmass = new((/ nRuns, nTracers, nTimes /), "float")
  dmass!0 = "run"
  dmass&run = runNames
  dmass!1 = "tracer"
  dmass&tracer = tracerNames
  dmass!2 = "time"
  dmass&time = history[:]&time
  dmass@units = "kg"
  
  ; Calculate the tracer masses for all the runs.
  do i = 0, nRuns - 1
    print("")
    print("")
    print("Run :  " + runNames(i))
    historyPath = basePath + runNames(i) + "/" + historyPattern + "/"
    test = calculateMass(historyPath, tracerNames)
    print(dimsizes(test))
    print(dimsizes(mass))
    print(i)
    mass(i, :, :) = (/test/)

    print("")
    do j = 0, nTracers - 1
      dmass(i, j, :) = mass(i, j, :) ;- mass(i, j, 0)
      
      print("  " + tracerNames(j) + " :  " + mass(i, j, nTimes-1) + " - " + \
        mass(i, j, 0) + " = " + dmass(i, j, nTimes-1))
    end do
  end do
  
  ; Create a panel plot of the mass difference.
  print("")
  print(" Plotting ...")
  pdffile="dmass_"+runNames(0)+"_"+historyPattern
  wks = gsn_open_wks("pdf", "$IMAGE_PATH/"+pdffile)
  plot = new(nRuns, graphic)
  
  res = True
;  res@gsnMaximize           = True
  res@vpHeightF             = 1.0 / nRuns              ; change aspect ratio
  res@vpWidthF              = 1.0
  res@gsnDraw = False
  res@gsnFrame = False
  
  res@trYMinF = 0
  res@trYMaxF = 5.0
  
;  res@tmYLFormat = "0*+^e"
;  res@tmYLPrecision = 3
  
  do i = 0, nRuns - 1
    res@gsnLeftString = runNames(i)
    res@gsnRightString = "Tg"
;    plot(i) = gsn_csm_xy(wks, dmass&time, dmass(i, nTracers - 1, :) / mass(i, nTracers-1, 0), res)
    plot(i) = gsn_csm_xy(wks, dmass&time/365., dmass(i, nTracers - 1, :) /1e9, res)
  end do
  
  resP = True
  resP@txString = "Total Mass"
  gsn_panel(wks, plot, (/ nRuns, 1 /), resP)
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

plotMass("/data4/mmills/","h0zm","b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg","SOOT01")
