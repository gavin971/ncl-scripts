load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

  filePath="/data5/mmills/b40.rcp4_5.2deg.wcm.carma.bc5tgIP.002/proc/"
  inFile="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.002.bcdep.nc"
  
  fileh=addfile(filePath+inFile,"r")
  
  
  bcphidry_test=fileh->bcphidry_test
  bcphiwet_test=fileh->bcphiwet_test
  bcphodry_test=fileh->bcphodry_test
  soot01dry_test=fileh->soot01dry_test
  soot01wet_test=fileh->soot01wet_test
  
  bcphidry_cntl=fileh->bcphidry_cntl
  bcphiwet_cntl=fileh->bcphiwet_cntl
  bcphodry_cntl=fileh->bcphodry_cntl
  
  time=fileh->time
  nTimes=dimsizes(time)
  
  data=new((/8,nTimes/),double)
  
  data(0,:)=cumsum(bcphidry_test,0)
  data(1,:)=cumsum(bcphiwet_test,0)
  data(2,:)=cumsum(bcphodry_test,0)
  
  data(3,:)=cumsum(bcphidry_cntl,0)
  data(4,:)=cumsum(bcphiwet_cntl,0)
  data(5,:)=cumsum(bcphodry_cntl,0)
  
  data(6,:)=cumsum(soot01dry_test,0)
  data(7,:)=cumsum(soot01wet_test,0)
  
  wks   = gsn_open_wks ("pdf","$IMAGE_PATH/BCdep5Tg") ; open workstation
  res                   = True                      ; plot mods desired
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperMargin         = 0.0          ; expand plot size to maximum
  
  res@gsnPaperOrientation = "landscape"
  res@xyLineColors     = (/"red","blue","black","red","blue","black","red","blue"/) ; change line color
  res@xyMonoLineColor  = False           ; different colors for different lines
  res@xyDashPatterns   = (/0,0,0,1,1,1,2,2/)
  res@xyLineThicknessF = 3.0 
  res@xyMonoLineThickness = True 
  
  res@trXMaxF=366.0d0
  res@trYMaxF=8.0d0
  res@trXMinF=1.0d0
  res@trYMinF=0.001d0
  
  res@trXLog=True
  res@trYLog=True
  
  res@tiXAxisString="Time (days)"
  res@tiYAxisString="Cumulative deposition (Tg)"
  
  
  
  plot = gsn_csm_xy (wks,time+1,data,res) ; draw plot
  
  delete(wks)
  
  mutter=systemfunc("mutt -a $IMAGE_PATH/BCdep5Tg.pdf mikejmills@mac.com < /dev/null")
  print(mutter)
  
  
