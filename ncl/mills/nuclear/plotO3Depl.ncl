load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("LOxLatMean")
function LOxLatMean(Rate[*][*]:float,j:integer)

local pi,wgt
begin
  pi=acos(0.0)*2.0
  if (j.eq.-1) then
    Rate!1 = "lat"
    lat=Rate&lat
    wgt=dble2flt(cos(lat*pi/180.))
    LOx=dim_avg_wgt_Wrap(Rate,wgt,0)
    delete(wgt)
  else
    LOx=Rate(:,j)
  end if
  delete(pi)
return(LOx)
end


undef("plotO3Depletion")
procedure plotO3Depletion(path:string, historyRun:string, controlRun:string)

begin
  print((/historyRun/))
  historyPattern = path+historyRun+"/proc/"+historyRun+".O3Col.h0zm.nc"
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
 
  controlPattern = path+controlRun+"/proc/"+controlRun+".O3Col.h0zm.nc"
  controlNames = systemfunc("ls " + controlPattern)
  control = addfiles(controlNames, "r")
  ListSetType(control, "cat")
  
  date_x=addfiles_GetVar(history, historyNames, "date")
  dx=dimsizes(date_x)
  print("dx:"+dx)
  print("date_x(0):"+date_x(0))
  print("date_x(dx-1):"+date_x(dx-1))
  date_c=addfiles_GetVar(control, controlNames, "date")
  dc=dimsizes(date_c)
  print("dc:"+dc)
  print("date_c(0):"+date_c(0))
  print("date_c(dc-1):"+date_c(dc-1))

  O3_col_x=addfiles_GetVar(history, historyNames, "O3_col")
  print("O3_col_x:"+dimsizes(O3_col_x))
  O3_col_c=addfiles_GetVar(control, controlNames, "O3_col")
  print("O3_col_c:"+dimsizes(O3_col_c))
  
  dc1=ind(date_c.eq.date_x(0))
  dc2=ind(date_c.eq.date_x(dx-1))
  
  print(dc1+": "+date_c(dc1))
  print(dc2+": "+date_c(dc2))
  
;  dc=dimsizes(O3_col_c)
;  dx=dimsizes(O3_col_x)
;  d1=min((/dc(0),dx(0)/))
 
;  O3_col_pct_loss=100.*(O3_col_x(0:d1-1,:)-O3_col_c)/O3_col_c
  O3_col_pct_loss=100.*(O3_col_x-O3_col_c(dc1:dc2,:))/O3_col_c(dc1:dc2,:)
  time=O3_col_x&time
;  print(time+" "+time/365.)
  time=time/365.  
  lat=O3_col_x&lat
  O3_col_pct_loss!1="lat"
  O3_col_pct_loss&lat=lat 
     
  wks = gsn_open_wks("pdf", "$IMAGE_PATH/O3Depletion."+historyRun)

  res = True
  resP = True
  resP@txString = "O3 Column % Change"
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
    
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.4           ; change aspect ratio of plot
  res@vpHeightF           = 0.55

;  res@gsnDraw = False
;  res@gsnFrame = False
  
;  res@trYMinF = 0
;  res@trYMaxF = 5.0
  
;  res@tmYLFormat = "0*+^e"
;  res@tmYLPrecision = 3
  res@pmLegendDisplayMode = "Always"
  nPanel=1
  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)
 

  res@tmXBMaxTicks     = 10
  res@tmXBMinorPerMajor= 11
;  res@tmXBMode         = "Manual"         ; Define own tick mark labels.
;  res@tmXBTickStartF     = 15/365.
;  res@tmXBTickSpacingF = 1
;  res@tmXBTickEndF     = 10
  res@tmXBAutoPrecision= False
  res@tmXBPrecision    = 1
  
;  res@tmXBValues         = major_ticks
;  res@tmXBLabels         = new(ndays,string)
;  res@tmXBLabels(:)     = ""

  
 ; Label mid tickmark with the month label.
;  res@tmXBLabels(::tk_per_mon)  = monthnames(months-1)

;  res@tmXBMinorOutwardLengthF = 0.02    ; Increase length and size of
;  res@tmXBMinorLengthF        = 0.01    ; minor tickmarks
;  res@tmXBMinorThicknessF     = 2.0
;  res@tmXBMajorThicknessF     = 2.0
;  res@tmXBMajorOutwardLengthF = 0.01
;  res@tmXBMajorLengthF        = 0.01
;===================================================================

  res@pmLegendSide           = "Top"               ; Change location of 
  res@pmLegendParallelPosF   = .4                  ; move units right
  res@pmLegendOrthogonalPosF = -0.4                ; move units down
  res@pmLegendWidthF         = 0.25                ; Change width and
  res@pmLegendHeightF        = 0.3                ; height of legend.
  res@lgPerimOn              = False               ; turn off box around
  res@lgLabelFontHeightF     = .04                 ; label font height
 
  res@tmYRLabelsOn           = True                ; Turn on right axis tick labels
  
  res@gsnLeftString = "O3 Column"
  res@gsnRightString = "% Change"
;  data= new((/5,dimsizes(time)/),float)
;  data(0,:)=LOxLatMean(O3_col_pct_loss,-1)
;  data(1,:)=O3_col_pct_loss(:,{86})
;  data(2,:)=O3_col_pct_loss(:,{66})
;  data(3,:)=O3_col_pct_loss(:,{46})
;  data(4,:)=O3_col_pct_loss(:,{26})
;  res@xyLineThicknessF = 2.0
;  res@xyLineColors      = (/"blue","red","green","black","orange"/); change line color
;  res@xyExplicitLegendLabels=(/"global mean", "86N", "66N", "46N", "26N"/)
;  plot(0) = gsn_csm_xy(wks, time, data, res)
  
;  data(0,:)=O3_col_pct_loss(:,{2})
;  data(1,:)=O3_col_pct_loss(:,{-26})
;  data(2,:)=O3_col_pct_loss(:,{-46})
;  data(3,:)=O3_col_pct_loss(:,{-66})
;  data(4,:)=O3_col_pct_loss(:,{-86})
;  delete(res@xyExplicitLegendLabels)
;  res@xyExplicitLegendLabels=(/"2N", "26S", "46S", "66S", "86S"/)
;  plot(1) = gsn_csm_xy(wks, time, data, res)
  
  print("create data")
  data= new((/4,dimsizes(time)/),float)
  print("created data")
;  data(0,:)=O3_col_pct_loss(:,{-86})
  data(0,:)=O3_col_pct_loss(:,{-46})
  data(1,:)=O3_col_pct_loss(:,{46})
  data(2,:)=O3_col_pct_loss(:,{86})
  print("calling LOxLatMean")
  data(3,:)=LOxLatMean(O3_col_pct_loss,-1)
  print("return from LOxLatMean")
  data!1="time"
  data&time=time
  res@xyLineThicknesses = (/4,4,4,8/)
  res@xyMonoLineThickness = False
  res@xyLineColors      = (/"blue","red","DarkOliveGreen","black"/); change line color
  res@xyDashPatterns=(/0,0,0,0,0/)
  res@xyExplicitLegendLabels=(/"46~S~o~N~S", "46~S~o~N~N", "86~S~o~N~N","global mean"/)
  res@trXMinF = time(0)
  res@trXMaxF = time(dimsizes(time)-1)
  res@trYMaxF = 10
  res@trYMinF = -70
  res@xyLineThicknessF = 3.0
;  res@tmXBLabelFontHeightF = 0.02
  print("plotting")
  plot(0) = gsn_csm_xy(wks, time, data, res)
  print("done plotting")
;  frame(wks)
  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(plot)

end

;path="$MODEL_PATH/ccsm/run/"
path="/data5/mmills/"
;historyRun="b40.rcp4_5.2deg.wcm.carma.bc5tgIndPak.001"
;controlRun="b40.rcp4_5.2deg.wcm.001"
historyRun="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.003"
controlRun="b40.rcp4_5.2deg.wcm.carma.bc0tg.002"

plotO3Depletion(path, historyRun, controlRun)
