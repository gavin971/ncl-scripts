load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;year_bounds=(/(/1960, 1979/),(/1979,2005/)/)

; WACCM4 1979-2005 ice extent
basepath="/glade/p/cesm/wawg/mmills/run/"
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

  ; Create a panel plot
  resP = True
  
  resP@txString = "Years after conflict"
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
  res@tiXAxisString = "Years after conflict"
  res@tiYAxisString = "Ice Extent"
    
  res2=res
  res2@xyLineThicknesses = (/0.5,0.5,0.5,0.5,3.0,3.0/)
  res2@xyMonoDashPattern = True
  res2@xyMonoLineColor   = False            ; same color for all lines
  res2@xyLineColors       = (/"lightblue","lightblue","pink","pink","midnightblue","red"/)
;  delete(res2@gsnXYFillColors)
  res2@gsnXYFillColors = (/"lightblue","transparent","pink","transparent","transparent"/)
  
  IMAGE_PATH=getenv("IMG_PATH")
  pdffile="iceExtentTimeSeries"
  wks   = gsn_open_wks ("pdf",IMAGE_PATH+"/"+pdffile) ; open workstation
  res@trYMinF = 8.5
  res@trYMaxF = 20
  plot1 = gsn_csm_xy (wks,time,dataNH,res2) ; draw plot

  res@trYMinF = 4
  res@trYMaxF = 24.5
  plot2 = gsn_csm_xy (wks,time,dataSH,res2) ; draw plot
  resP@tiYAxisFontHeightF=0.015
  newplot = gsn_attach_plots(plot1,plot2,res,resP)
  draw(plot1)
  
delete(wks)
  
