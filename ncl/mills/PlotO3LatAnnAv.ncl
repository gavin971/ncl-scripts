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


undef("plotO3LatAnnAv")
procedure plotO3LatAnnAv(historyRun:string, controlRun:string)

begin
  print((/historyRun/))
  historyPattern = historyRun+"/"+historyRun+".mean.h3.1996.nc"
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
 
  controlPattern = controlRun+"/"+controlRun+".mean.h3.1996.nc"
  controlNames = systemfunc("ls " + controlPattern)
  control = addfiles(controlNames, "r")
  ListSetType(control, "cat")

  O3_col_x=addfiles_GetVar(history, historyNames, "O3_col")
  O3_col_c=addfiles_GetVar(control, controlNames, "O3_col")
  print("O3_col_x:"+dimsizes(O3_col_x))
  print("O3_col_c:"+dimsizes(O3_col_c))
  dc=dimsizes(O3_col_c)
  dx=dimsizes(O3_col_x)
  d1=min((/dc(0),dx(0)/))
  
  O3_col_c=O3_col_c/2.687e16 ; Convert to Dobson Units
  O3_col_c@units="Dobson Units"
  O3_col_x=O3_col_x/2.687e16 ; Convert to Dobson Units
  O3_col_x@units="Dobson Units" 
  
  time=O3_col_c&time
  lat=O3_col_x&lat

  day1=31+28+31+30+15  ; May 15
  
  days=day1+time
  
;  jan1=(/366-day1,366-day1+365,366-day1+2*365,366-day1+3*365/)
       
;  O3ColXAnnAv=new((/dimsizes(lat)/),float)
;  O3ColCAnnAv=new((/dimsizes(lat)/),float)
  
  O3ColXAnnAv=dim_avg_Wrap(O3_col_x(lat|:,time|:))
  O3ColCAnnAv=dim_avg_Wrap(O3_col_c(lat|:,time|:))
  
  O3ColXAnnMin=O3ColXAnnAv ; get metadata
  O3ColCAnnMin=O3ColCAnnAv ; get metadata
  O3ColXAnnMax=O3ColXAnnAv ; get metadata
  O3ColCAnnMax=O3ColCAnnAv ; get metadata
  
  O3ColXAnnMin=dim_min(O3_col_x(lat|:,time|:))
  O3ColCAnnMin=dim_min(O3_col_c(lat|:,time|:))

  O3ColXAnnMax=dim_max(O3_col_x(lat|:,time|:))
  O3ColCAnnMax=dim_max(O3_col_c(lat|:,time|:))
  
  wks = gsn_open_wks("pdf", "O3LatAnnAv.h3")

  res = True
  resP = True
;  resP@txString = "O3 Column % Change"
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
    
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.9

;  res@pmLegendDisplayMode = "Always"
;  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
;  plot = new(nPanel, graphic)
 

;  res@tmYLMode           = "Manual" 
;  res@tmYLTickStartF     = -70
;  res@tmYLTickSpacingF        = 10
;  res@tmYRMode           = "Manual"     
;  res@tmYRTickStartF     = -70
;  res@tmYRTickSpacingF   = 10
;  res@gsnCenterString      = "Annual average ozone column (Dobson Units)"
  res@tiYAxisString      = "Ozone column (DU)"
  res@tiXAxisString      = "Latitude"
  
  res@pmLegendSide           = "Top"               ; Change location of 
  res@pmLegendParallelPosF   = .45                  ; move units right
  res@pmLegendOrthogonalPosF = -0.5                ; move units down
;  res@lgLabelAngleF          = 290
;  res@lgOrientation = "horizontal"                 ; swap the orientation of the legend
  res@pmLegendWidthF         = 0.25                ; Change width and
  res@pmLegendHeightF        = 0.15                ; height of legend.
  res@lgPerimOn              = False               ; turn off box around
  res@lgLabelFontHeightF     = .025                 ; label font height
 
  res@tmYRLabelsOn           = True                ; Turn on right axis tick labels
    
  res@xyMonoLineThickness = False
  res@xyLineThicknesses = (/4.0,4.0,2.0,2.0,2.0,2.0/)
  res@xyLineColors      = (/"red","blue","brown","dodgerblue","brown","dodgerblue"/); change line color
  res@xyDashPatterns=(/0,0,2,2,2,2/)
  res@xyExplicitLegendLabels=(/"soot", "control"/)
;  res@trXMaxF = time(dimsizes(time)-1)
;  res@trYMaxF = 15
  res@trYMinF = 0
  res@tmXBLabelFontHeightF = 0.03
  res@tmYLLabelFontHeightF = 0.03

  data= new((/6,dimsizes(lat)/),float)
  data(0,:)=O3ColXAnnAv(:)
  data(1,:)=O3ColCAnnAv(:)
  data(2,:)=O3ColXAnnMin(:)
  data(3,:)=O3ColCAnnMin(:)
  data(4,:)=O3ColXAnnMax(:)
  data(5,:)=O3ColCAnnMax(:)
 plot = gsn_csm_xy(wks, lat, data, res)

;  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
;  delete(plot)

end

plotO3LatAnnAv("geo2","nongeo2")
