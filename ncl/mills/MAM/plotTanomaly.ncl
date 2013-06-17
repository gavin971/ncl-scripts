load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;runNames = (/ "b.e11.B55TRW5CN.f19_g16.Pin16-18.dgnumhi1um.001", \
;              "b.e11.B55TRW5CN.f19_g16.Pin20-30.dgnh1um.noqbo.001" , \
;              "b.e11.B55TRW5CN.f19_g16.Pin17-25.dgnh1um.noqbo.001" /)
runNames = (/ "b.e11.B55TRW5CN.f19_g16.Pin16-20.dgnh400nm.noqbo.001", \
              "b.e11.B55TRW5CN.f19_g16.Pin16-20.dgnh1um.noqbo.001" , \
              "b.e11.B55TRW5CN.f19_g16.Pin16-20.dgnh10um.noqbo.001" /)
nRuns = dimsizes(runNames)

cntlRun  =    "b.e11.B55TRW5CN.f19_g16.MJMtest.013"          
;cntlRunNames  =   (/"b.e11.B55TRW5CN.f19_g16.novolc.dgnh400nm.noqbo.001" , \
;              "b.e11.B55TRW5CN.f19_g16.MJMtest.013c" , \
;              "b.e11.B55TRW5CN.f19_g16.novolc.dgnh10um.noqbo.001" /)         
;rootpath =    "/glade/scratch/mmills/archive/"
rootpath="/glade/p/cesm/wawg/mmills/run/"

levels = (/ 30., 50., 70., 100., 1013.25 /)
nLevs = dimsizes(levels)

filename=cntlRun+".cam.h0.tropicsavg.nc"
filepath=rootpath+cntlRun+"/atm/proc/h0tm/"

fileid=addfile(filepath+filename,"r")
date=fileid->date
d1=ind(date.eq.19910201)
d2=ind(date.eq.19930101)
nTimes=d2-d1+1
tempcntl=fileid->T(d1:d2,:)
TScntl=fileid->TS(d1:d2)
datecntl=date(d1:d2)
time=fileid->time(d1:d2)
;print(time)
time=time-fileid->time(d1-1)
time@units="days since "+date(d1-1)+" 00:00:00"
;print(time)
delete(date)
delete(fileid)

temp = new((/ nLevs, nRuns+1, nTimes /), "float")
temp!0 = "lev"
temp&lev = levels
temp!1 = "run"
temp!2 = "time"
temp&time=time-15      
temp@units = "K"
temp(:,:,:)=0.0e0
print(temp&time)

obsfile="/glade/p/cesm/wawg/mmills/inputdata/rich_gridded_2009.nc"
fileid=addfile(obsfile,"r")
t1=floattoint(365.25*(1991-1958)+0.5)
t2=floattoint(365.25*(1993-1958)-30.5)
anomalies=fileid->anomalies ; ( time, pressure, lat, lon )
anomalies@_FillValue = -1e+30
tobs=dim_avg_Wrap(anomalies(:,:,:,:))
tobsTropics=dim_avg_Wrap(tobs({t1:t2},:,{-25:25}))
print("   Max obs:"+max(tobsTropics))
print("   Min obs:"+min(tobsTropics))
; Normalize against 1989-1990
t3=floattoint(365.25*(1989-1958)+0.5)
t4=floattoint(365.25*(1991-1958)-30.5)
tobsNormTime=dim_avg_Wrap(tobs({t3:t4},:,{-25:25}))
tobsNorm=dim_avg_n_Wrap(tobsNormTime,0)
;tobsTropics=tobsTropics-tobsNorm
print("   Max obs:"+max(tobsTropics))
print("   Min obs:"+min(tobsTropics))

pdffile="Tanomaly_Pinatubo"
wks = gsn_open_wks("pdf", "$IMAGE_PATH/"+pdffile)

do l = 0, nLevs - 1
  print("")
  print("Level :  " + levels(l))
  do i = 0, nRuns - 1
    print("")
    print("  Run :  " + runNames(i))
    filename=runNames(i)+".cam.h0.tropicsavg.nc"
    filepath=rootpath+runNames(i)+"/atm/proc/h0tm/"

    fileid=addfile(filepath+filename,"r")
    date=fileid->date
    nt=dimsizes(date)
;    print("  nt="+nt)
    d1=ind(date.eq.19910201)
    d2=ind(date.eq.19930101)

    if (levels(l).eq.1013.25) then
      temp(l, i, :) = (/fileid->TS(d1:d2)/)
      temp(l, i, :) = temp(l, i, :) - (/TScntl(:)/)
    else
      temp(l, i, :) = (/fileid->T(d1:d2,{levels(l)})/)
      temp(l, i, :) = temp(l, i, :) - (/tempcntl(:,{levels(l)})/)
    end if
    
    delete(date)
  end do
  if (levels(l).ne.1013.25) then
    temp(l, nRuns, :) = (/tobsTropics(:,{levels(l)})/) - tobsNorm({levels(l)})
  end if
;  print("   Max obs:"+max(temp(l, nRuns, :)))
;  print("   Min obs:"+min(temp(l, nRuns, :)))
  res = True
  res@gsnMaximize           = True
  res@vpHeightF             = 0.5              ; change aspect ratio
  res@vpWidthF              = 1.0
  res@xyLineColors    = (/"blue","green4","red","black"/)
  res@xyMonoLineColor = False
  res@xyMonoDashPattern = True
  res@xyMonoLineThickness = False
  res@xyLineThicknesses = (/4.0, 4.0, 4.0, 8.0/)
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
  res@trYMinF = -2.0
  if (l .ge. 3) then
    delete(res@trYMinF)
  end if
  res@tmXMajorGrid = True
  res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  res@tmYMajorGrid = True
  res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
;  print("   Max all runs:"+max(temp(l, :, :)))
;  print("   Min all runs:"+min(temp(l, :, :)))
  res@gsnLeftString = " "
  res@gsnCenterString = levels(l)+" hPa"
  res@gsnRightString = " "
  res@tiXAxisString = "Year"
  res@tiYAxisString = "Temperature anomaly (K)"
;  print(temp&time)
  plot = gsn_csm_xy(wks, 1991+temp&time/365., temp(l, :, :), res)
  
end do

delete(wks)


