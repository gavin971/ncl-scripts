load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;runNames = (/ "b.e11.B55TRW5CN.f19_g16.Pin16-18.dgnumhi1um.001", \
;              "b.e11.B55TRW5CN.f19_g16.Pin20-30.dgnh1um.noqbo.001" , \
;              "b.e11.B55TRW5CN.f19_g16.Pin17-25.dgnh1um.noqbo.001" /)
;runNames = (/ "b.e11.B55TRW5CN.f19_g16.Pin16-20.dgnh400nm.noqbo.001", \
;              "b.e11.B55TRW5CN.f19_g16.Pin16-20.dgnh1um.noqbo.001" , \
;              "b.e11.B55TRW5CN.f19_g16.Pin16-20.dgnh10um.noqbo.001" /)
;runNames = (/"b.e12.B55TRW5CN.f19_g16.course.001", \
;             "b.e12.B55TRW5CN.f19_g16.CoursePin.002" /)
;runNames = (/"b.e12.B55TRW5CN.f19_g16.CoursePin16-18.001", \
;             "b.e12.B55TRW5CN.f19_g16.CoursePin16-18.002" /)
;runNames = (/ "b.e12.B55TRW5CN.f19_g16.CoursePin17-25.QBOon.001", \
;            "b.e12.B55TRW5CN.f19_g16.CoursePin17-25.QBOon.002" /)
runNames = (/ "b.e12.B55TRW5CN.f19_g16.Pin10Tg.16-20km.QBOon.001", \
              "b.e12.B55TRW5CN.f19_g16.Pin10Tg.16-20km.QBOon.002" /)
nRuns = dimsizes(runNames)

;cntlRun  =    "b.e11.B55TRW5CN.f19_g16.MJMtest.013"          
;cntlRunNames  =   (/"b.e11.B55TRW5CN.f19_g16.novolc.dgnh400nm.noqbo.001" , \
;              "b.e11.B55TRW5CN.f19_g16.MJMtest.013c" , \
;              "b.e11.B55TRW5CN.f19_g16.novolc.dgnh10um.noqbo.001" /)         
;cntlRunNames  =    (/ "b.e12.B55TRW5CN.f19_g16.CourseBkgd.001" , \
;                      "b.e12.B55TRW5CN.f19_g16.CourseBkgd.002" /)
cntlRunNames  =    (/ "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.001" , \
                      "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.002" /)
rootpath =    "/glade/scratch/mmills/archive/"
;rootpath="/glade/p/cesm/wawg/mmills/run/"

year1=1991
year2=1994
;year2=1993
date1=year1*10000+0201
date2=year2*10000+0701

levels = (/ 30., 50., 70., 100., 1013.25 /)
nLevs = dimsizes(levels)

cntlRun = cntlRunNames(0)

filename=cntlRun+".cam.h0.tropicsavg.nc"
filepath=rootpath+cntlRun+"/atm/proc/h0tm/"

fileid=addfile(filepath+filename,"r")
date=fileid->date
d1=ind(date.eq.date1)
d2=ind(date.eq.date2)
nTimes=d2-d1+1
;tempcntl=fileid->T(d1:d2,:)
;TScntl=fileid->TS(d1:d2)
datecntl=date(d1:d2)
time=fileid->time(d1:d2)
;print(time)
time=time-fileid->time(0)
time@units="days since "+date(0)+" 00:00:00"
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
t1=floattoint(12*(year1-1958))
t2=floattoint(12*(year2-1958)-1)
t2=t2+6
anomalies=fileid->anomalies ; ( time, pressure, lat, lon )
anomalies@_FillValue = -1e+30
latobs = anomalies&lat
gwobs = latRegWgt(latobs,"float",0) 
gwobs&lat=latobs
tobsTropicsAll=wgt_areaave_Wrap(anomalies(:,:,{-25:25},:),gwobs({-25:25}),1.0,0)
;print("   Max obs:"+max(tobsTropics))
;print("   Min obs:"+min(tobsTropics))
; Normalize against 1986-1990
t3=floattoint(12*(1986-1958))
t4=floattoint(12*(1990-1958))
tobsNorm=new((/12,12/),"float")
printVarSummary(tobsTropicsAll)
print("t3="+t3+", t4="+t4)
do mon=0,11 
  do lev=0,11
    tobsNorm(mon,lev)=avg(tobsTropicsAll(t3+mon:t4+mon:12,lev))
  end do
end do
tobsTropics=tobsTropicsAll(t1:t2,:)
do yr = 0,year2-1-year1
  tobsTropics(yr*12:yr*12+11,:) = tobsTropics(yr*12:yr*12+11,:) - tobsNorm(:,:)
end do
delete(fileid)
delete(latobs)
delete(gwobs)
print("   Max obs:"+max(tobsTropics))
print("   Min obs:"+min(tobsTropics))

surfObsFile="/glade/p/cesm/wawg/mmills/inputdata/HadCRUT3v.nc"
fileid=addfile(surfObsFile,"r")
var_in = fileid->temp(:,0,:,:)
delete(fileid)
latobs = var_in&latitude
gwobs = latRegWgt(latobs,"float",0) 
var_w = wgt_areaave_Wrap(var_in,gwobs,1.0,0)
t1=(year1-1850)*12
t2=(year2-1850)*12-1
t2=t2+6
TSobsGlobal = var_w(t1:t2)
; Normalize against 1986-1990
t3=(1986-1850)*12
t4=(1990-1850)*12
TSnorm=new(12,"float")
do mon=0,11 
  TSnorm(mon)=avg(var_w(t3+mon:t4+mon:12))
end do
do yr = 0,year2-1-year1
  TSobsGlobal(yr*12:yr*12+11) = TSobsGlobal(yr*12:yr*12+11) - TSnorm(0:11)
end do
delete(var_w)

pdffile="Tanomaly_Pinatubo"
pdfpath="$IMAGE_PATH/"+runNames(0)+"/"
mkdir=systemfunc("mkdir -p "+pdfpath)
wks = gsn_open_wks("pdf", pdfpath+pdffile)

do l = 0, nLevs - 1
  print("")
  print("Level :  " + levels(l))
  do i = 1, nRuns
    print("")
    print("  Expt Run :  " + runNames(i-1))
    print("  Cntl Run :  " + cntlRunNames(i-1))
    
    if (levels(l).eq.1013.25) then

      filename=runNames(i-1)+".cam.h0.globavg.nc"
      filepath=rootpath+runNames(i-1)+"/atm/proc/h0gm/"
      fileid=addfile(filepath+filename,"r")

      cntlFilename=cntlRunNames(i-1)+".cam.h0.globavg.nc"
      cntlFilepath=rootpath+cntlRunNames(i-1)+"/atm/proc/h0gm/"
      cntlFileid=addfile(cntlFilepath+cntlFilename,"r")
    
      T      = fileid->TS
      Tcntl = cntlFileid->TS
      
;      temp(l, i, :) = (/fileid->TS(d1:d2)/)
;      TScntl   = cntlFileid->TS(d1:d2)
;      temp(l, i, :) = temp(l, i, :) - (/TScntl(:)/)
    else

      filename=runNames(i-1)+".cam.h0.tropicsavg.nc"
      filepath=rootpath+runNames(i-1)+"/atm/proc/h0tm/"
      fileid=addfile(filepath+filename,"r")

      cntlFilename=cntlRunNames(i-1)+".cam.h0.tropicsavg.nc"
      cntlFilepath=rootpath+cntlRunNames(i-1)+"/atm/proc/h0tm/"
      cntlFileid=addfile(cntlFilepath+cntlFilename,"r")
      
      T      = fileid->T(:,{levels(l)})
      Tcntl = cntlFileid->T(:,{levels(l)})
;      temp(l, i, :) = (/fileid->T(d1:d2,{levels(l)})/)
;      tempcntl = cntlFileid->T(d1:d2,{levels(l)})
;      temp(l, i, :) = temp(l, i, :) - (/tempcntl(:)/)
    end if
    
    date=fileid->date
    d1=ind(date.eq.date1)
    d2=ind(date.eq.date2)
    
    cntldate=fileid->date
    cd1=ind(cntldate.eq.date1)
    cd2=ind(cntldate.eq.date2)
    
    temp(l, i, :) = (/T(d1:d2)/) - (/Tcntl(cd1:cd2)/)
    delete(T)
    delete(Tcntl)
    delete(cntldate)
    
    delete(date)
  end do
  if (levels(l).eq.1013.25) then
;    printVarSummary(temp)
;    printVarSummary(TSobsGlobal)
;    printVarSummary(TSnorm)
    temp(l, 0, :) = (/TSobsGlobal(:)/)
  else
    temp(l, 0, :) = (/tobsTropics(:,{levels(l)})/)
  end if
;  print("   Max obs:"+max(temp(l, nRuns, :)))
;  print("   Min obs:"+min(temp(l, nRuns, :)))
  res = True
  res@gsnMaximize           = True
  res@vpHeightF             = 0.5              ; change aspect ratio
  res@vpWidthF              = 1.0
  res@xyLineColors    = (/"black","blue","green4","red"/)
  res@xyMonoLineColor = False
  res@xyMonoDashPattern = True
  res@xyMonoLineThickness = False
  res@xyLineThicknesses = (/8.0, 4.0, 4.0, 4.0/)
  res@trXMinF = year1
  res@trXMaxF = year2
  res@tmXBMinorOn = False
;  res@tmXBMinorPerMajor = 11
  nticks=(year2-year1)*12+1
  res@tmXBMaxTicks = nticks
  res@tmXBMode= "Explicit"
  res@tmXBValues =fspan(year1,year2,nticks)
  res@tmXBLabelFontHeightF = 0.015
;  XBLabels = new(nticks,"string")
  res@tmXBLabels = (/"91","F","M","A","M","J","J","A","S","O","N","D", \
                     "92","F","M","A","M","J","J","A","S","O","N","D", \
                     "93","F","M","A","M","J","J","A","S","O","N","D", \
                     "94" /)
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
  if (levels(l).eq.1013.25) then
    res@gsnCenterString = "Surface, seasonally adjusted global average"
  else
    res@gsnCenterString = levels(l)+" hPa, seasonally adjusted tropical average"
  end if
  res@gsnRightString = " "
  res@tiXAxisString = "Year"
  res@tiYAxisString = "Temperature anomaly (K)"
;  print(temp&time)
  plot = gsn_csm_xy(wks, 1991+temp&time/365., temp(l, :, :), res)
  
end do

delete(wks)


