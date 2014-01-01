; This script creates Figure 4 for 2013 paper "Multi-decadal global cooling 
; and unprecedented ozone loss following a regional nuclear conflict" by 
; Mills, Toon, Lee-Taylor, and Robock. 

; Figure 4 caption: Change in sea ice extent (%) for the 5 Tg BC experiment 
; ensemble range is shown relative to the control ensemble average. Sea ice 
; extent is defined as the area of all sea surface gridpoints with ice fraction 
; greater than 15%. Red is for the Southern Hemisphere average, and blue is for 
; the Northern Hemisphere average.

; Script written by Mike Mills
; Modified September 24, 2013

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;year_bounds=(/(/1960, 1979/),(/1979,2005/)/)

ostype = systemfunc("echo $OSTYPE")
if (ostype.eq."darwin") then
  basepath="$MODEL_PATH/cesm/runs/cesm1_0/"
else
  basepath="/glade/p/cesm/wawg/mmills/run/"
end if

; WACCM ice extent
casesExpt=(/"b40.rcp4_5.2deg.wcm.carma.bc5tgIP.003", \
            "b40.rcp4_5.2deg.wcm.carma.bc5tgIP.006", \
            "b40.rcp4_5.2deg.wcm.carma.bc5tgIP.007"/)
casesCntl=(/"b40.rcp4_5.2deg.wcm.carma.bc0tg.006", \
            "b40.rcp4_5.2deg.wcm.carma.bc0tg.007", \
            "b40.rcp4_5.2deg.wcm.carma.bc0tg.008"/)
            
ncases=dimsizes(casesExpt)
delete(filesExpt)
filesExpt=basepath+casesExpt+"/ice/proc/iextent_reg_avg_timeseries.dat"
delete(filesCntl)
filesCntl=basepath+casesCntl+"/ice/proc/iextent_reg_avg_timeseries.dat"
ncol=3
ignorelines=1

do i = 0, ncases-1
  print("reading "+filesExpt(i))
  iextent_table=readAsciiTable(filesExpt(i),ncol,"float",ignorelines)
  if (i.eq.0) then
    ndims=dimsizes(iextent_table)
    date=iextent_table(:,0)
    iextNHExpt=new(ndims,float)
    iextNHCntl=new(ndims,float)
    iextSHExpt=new(ndims,float)
    iextSHCntl=new(ndims,float)
  end if 
  iextNHExpt(:,i)=iextent_table(:,2)/1e12
  iextSHExpt(:,i)=iextent_table(:,1)/1e12
  delete(iextent_table)
  print("reading "+filesCntl(i))
  iextent_table=readAsciiTable(filesCntl(i),ncol,"float",ignorelines)
  iextNHCntl(:,i)=iextent_table(:,2)/1e12
  iextSHCntl(:,i)=iextent_table(:,1)/1e12
  delete(iextent_table)  
end do

time=date-date(0)

dataNH=new((/6,ndims(0)/),float)
dataSH=new((/6,ndims(0)/),float)
dataRel=new((/6,ndims(0)/),float)

dataNH(0,:)=dim_min(iextNHExpt)
dataNH(1,:)=dim_max(iextNHExpt)
dataNH(2,:)=dim_min(iextNHCntl)
dataNH(3,:)=dim_max(iextNHCntl)
dataNH(4,:)=dim_avg(iextNHExpt)
dataNH(5,:)=dim_avg(iextNHCntl)

dataSH(0,:)=dim_min(iextSHExpt)
dataSH(1,:)=dim_max(iextSHExpt)
dataSH(2,:)=dim_min(iextSHCntl)
dataSH(3,:)=dim_max(iextSHCntl)
dataSH(4,:)=dim_avg(iextSHExpt)
dataSH(5,:)=dim_avg(iextSHCntl)

dataRel(0,:)=dataNH(0,:)/dataNH(5,:)
dataRel(1,:)=dataNH(1,:)/dataNH(5,:)
dataRel(2,:)=dataSH(0,:)/dataSH(5,:)
dataRel(3,:)=dataSH(1,:)/dataSH(5,:)
dataRel(4,:)=dataNH(4,:)/dataNH(5,:)
dataRel(5,:)=dataSH(4,:)/dataSH(5,:)

dataRel=(dataRel-1.)*100.

  ; Create a panel plot
  resP = True
  
  resP@txString = "Time (years since BC injection)"
  resP@gsnAttachPlotsXAxis =True
  resP@gsnPaperOrientation = "landscape"
  resP@gsnMaximize         = True          ; expand plot size to maximum
  resP@tmYMajorGrid = True
 
  
  res                   = True                      ; plot mods desired
  res@tmYMajorGrid = True
  res@tmYMajorGridLineDashPattern = 2
  res@vpWidthF            = 1.3           ; change aspect ratio of plot
  res@vpHeightF           = 0.3
  res@gsnDraw = False
  res@gsnFrame     = False                        ; don't advance frame yet
  res@xyLineThicknessF = 5.0
  
  res@trXMinF = 0.0
  res@trXMaxF = max(time)
  res@tiXAxisString = "Time (years)"
  res@tiYAxisString = "Ice Extent"
    
  res2=res
  res2@xyLineThicknesses = (/0.5,0.5,0.5,0.5,3.0,3.0/)
  res2@xyMonoDashPattern = True
  res2@xyMonoLineColor   = False            ; same color for all lines
  res2@xyLineColors       = (/"lightblue","lightblue","pink","pink","midnightblue","red"/)
;  delete(res2@gsnXYFillColors)
  res2@gsnXYFillColors = (/"lightblue","transparent","pink","transparent","transparent"/)
  
  IMAGE_PATH=getenv("IMG_PATH")
  pdffile="iceExtentTimeSeriesBC5Tg"
  wks   = gsn_open_wks ("pdf",IMAGE_PATH+"/"+pdffile) ; open workstation
  res@trYMinF = 8.5
  res@trYMaxF = 20.0
  plot1 = gsn_csm_xy (wks,time,dataNH,res2) ; draw plot

  res@trYMinF = 4.0
  res@trYMaxF = 24.5
  plot2 = gsn_csm_xy (wks,time,dataSH,res2) ; draw plot
  resP@tiYAxisFontHeightF=0.015
  newplot = gsn_attach_plots(plot1,plot2,res,resP)
  draw(plot1)
  
  delete(wks)
  pdffile="Fig05_iceExtentTimeSeriesPercentBC5Tg"
  wks   = gsn_open_wks ("pdf",IMAGE_PATH+"/nuclear/"+pdffile) ; open workstation  

;  frame(wks)
  
  delete(res)
  res                   = True
  res@gsnDraw      = False                    ; Turn off draw and frame so
  res@gsnFrame     = False                    ; we can attach some text.
  res@tmYMajorGrid = True
  res@tmYMajorGridLineDashPattern = 2
  res@vpWidthF            = 0.7           ; change aspect ratio of plot
  res@vpHeightF           = 0.5
  res@xyLineThicknessF = 5.0
  res@trXMinF = 0.0
  res@trXMaxF = max(time)
  res@trYMinF = -15
  res@trYMaxF = 85
  res@tiXAxisString = "Time (years since BC injection)"
  res@tiYAxisString = "Ice Extent Change (%)"
  res@xyLineThicknesses = (/0.5,0.5,0.5,0.5,3.0,3.0/)
  res@xyMonoDashPattern = True
  res@xyMonoLineColor   = False            ; same color for all lines
  res@xyLineColors      = (/"lightblue","lightblue","pink","pink","midnightblue","red"/)
  res@gsnXYFillColors   = (/"lightblue","transparent","pink","transparent","transparent"/)
  
  plot = gsn_csm_xy (wks,time,dataRel,res)
  
;***********************************************
; legend resources
;***********************************************
  lgres                    = True
  lgres@lgLineColors     = (/"midnightblue","red"/)
  lgres@lgMonoLineLabelFontColor = False
  lgres@lgLineLabelFontColors = (/"midnightblue","red"/)
  lgres@lgMonoDashIndex  = True
  lgres@lgLineThicknessF   = res@xyLineThicknessF  ; legend line thickness = contour line thickness
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .2            ; set the legend label font thickness
  lgres@vpWidthF           = 0.15           ; width of legend (NDC)
  lgres@vpHeightF          = 0.15            ; height of legend (NDC)
  lgres@lgPerimOn = False
;  lgres@lgPerimColor       = "orange"       ; draw the box perimeter in orange
;  lgres@lgPerimThicknessF  = 5.0            ; thicken the box perimeter
  lbid = gsn_create_legend(wks,2,(/" NH"," SH"/),lgres)         ; create legend

  amres = True
  amres@amParallelPosF   = 0.30		    ; move legend to the right
  amres@amOrthogonalPosF = -0.25             ; move the legend up
  annoid1 = gsn_add_annotation(plot,lbid,amres)   ; add legend to plot

  draw(plot)
  frame(wks)
  
delete(wks)
  
