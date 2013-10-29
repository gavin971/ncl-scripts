load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

host   = systemfunc("hostname")
ostype = systemfunc("echo $OSTYPE")

if (isStrSubset(host,"caldera").or.isStrSubset(host,"geyser").or.isStrSubset(host,"yslogin")) then
  rootpath =    "/glade/scratch/mmills/archive/"
  ;rootpath="/glade/p/cesm/wawg/mmills/run/"
else
  rootpath =    "$MODEL_PATH/cesm/case/"
end if

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

;tag="15Tg.16-20km.1col"
;runNames = (/ "b.e12.B55TRW5CN.f19_g16.Pin15Tg.16-20km.QBOon.001", \
;              "b.e12.B55TRW5CN.f19_g16.Pin15Tg.16-20km.QBOon.002", \
;              "b.e12.B55TRW5CN.f19_g16.Pin15Tg.16-20km.QBOon.003" /)

;tag="10Tg.16-20km.1col"
;runNames = (/ "b.e12.B55TRW5CN.f19_g16.Pin10Tg.16-20km.QBOon.001", \
;              "b.e12.B55TRW5CN.f19_g16.Pin10Tg.16-20km.QBOon.002", \
;              "b.e12.B55TRW5CN.f19_g16.Pin10Tg.16-20km.QBOon.003" /)

;runNames = (/ "b.e12.B55TRW5CN.f19_g16.Pin10Tg.18-22km.QBOon.001", \
;              "b.e12.B55TRW5CN.f19_g16.Pin10Tg.18-22km.QBOon.002", \
;              "b.e12.B55TRW5CN.f19_g16.Pin10Tg.18-22km.QBOon.003" /)

;runNames = (/ "b.e12.B55TRW5CN.f19_g16.Pin15Tg.20-27km.0-18N.96-118E.QBOon.001", \
;              "b.e12.B55TRW5CN.f19_g16.Pin15Tg.20-27km.0-18N.96-118E.QBOon.002", \
;              "b.e12.B55TRW5CN.f19_g16.Pin15Tg.20-27km.0-18N.96-118E.QBOon.003" /)

tag="20Tg.20-27km.100col"
;runNames = (/ \
;   "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-27km.0-18N.96-118E.QBOon.001", \
;   "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-27km.0-18N.96-118E.QBOon.002", \
;   "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-27km.0-18N.96-118E.QBOon.003", \
;   "b.e12.B55TRW5CN.f19_g16.Pin15Tg.20-27km.0-18N.96-118E.QBOon.001", \
;   "b.e12.B55TRW5CN.f19_g16.Pin15Tg.20-27km.0-18N.96-118E.QBOon.002", \
;   "b.e12.B55TRW5CN.f19_g16.Pin15Tg.20-27km.0-18N.96-118E.QBOon.003", \
;   "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-27km.0-18N.96-118E.QBOon.001", \
;   "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-27km.0-18N.96-118E.QBOon.002", \
;   "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-27km.0-18N.96-118E.QBOon.003", \
;   "b.e12.B55TRW5CN.f19_g16.Pin10Tg.16-20km.QBOon.001", \
;   "b.e12.B55TRW5CN.f19_g16.Pin10Tg.16-20km.QBOon.002", \
;   "b.e12.B55TRW5CN.f19_g16.Pin10Tg.16-20km.QBOon.003" \
;               /)
runNames = (/ \
   "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-27km.0-18N.96-118E.QBOon.avg", \
   "b.e12.B55TRW5CN.f19_g16.Pin15Tg.20-27km.0-18N.96-118E.QBOon.avg", \
   "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-27km.0-18N.96-118E.QBOon.avg", \
   "b.e12.B55TRW5CN.f19_g16.Pin10Tg.16-20km.QBOon.avg", \
   "b.e12.B55TRW5CN.f19_g16.Pin20Tg.16-18km.1col.QBOon.avg" \
               /)

;tag="20Tg.20-27km.100col"
;runNames = (/ \
;   "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-27km.0-18N.96-118E.QBOon.001", \
;   "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-27km.0-18N.96-118E.QBOon.002", \
;   "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-27km.0-18N.96-118E.QBOon.003", \
;   "b.e12.B55TRW5CN.f19_g16.Pin15Tg.20-27km.0-18N.96-118E.QBOon.001", \
;   "b.e12.B55TRW5CN.f19_g16.Pin15Tg.20-27km.0-18N.96-118E.QBOon.002", \
;   "b.e12.B55TRW5CN.f19_g16.Pin15Tg.20-27km.0-18N.96-118E.QBOon.003", \
;   "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-27km.0-18N.96-118E.QBOon.001", \
;   "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-27km.0-18N.96-118E.QBOon.002", \
;   "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-27km.0-18N.96-118E.QBOon.003" \
;               /)

;tag="10Tg.20-27km.100col"
;runNames = (/ "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-27km.0-18N.96-118E.QBOon.001", \
;              "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-27km.0-18N.96-118E.QBOon.003" /)

nRuns = dimsizes(runNames)

;cntlRun  =    "b.e11.B55TRW5CN.f19_g16.MJMtest.013"          
;cntlRunNames  =   (/"b.e11.B55TRW5CN.f19_g16.novolc.dgnh400nm.noqbo.001" , \
;              "b.e11.B55TRW5CN.f19_g16.MJMtest.013c" , \
;              "b.e11.B55TRW5CN.f19_g16.novolc.dgnh10um.noqbo.001" /)         
;cntlRunNames  =    (/ "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.001" , \
;                      "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.003" /)
;cntlRunNames  =    (/ "b.e12.B55TRW5CN.f19_g16.CourseBkgd.001" , \
;                      "b.e12.B55TRW5CN.f19_g16.CourseBkgd.002" /)
cntlRunNames  =    (/ \
   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.avg", \
   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.avg", \
   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.avg", \
   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.avg", \
   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.avg" \
   /)
;   , \
;   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.001", \
;   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.002", \
;   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.003", \
;   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.001", \
;   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.002", \
;   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.003", \
;   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.001", \
;   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.002", \
;   "b.e12.B55TRW5CN.f19_g16.CourseBkgd.QBOon.003" \
;                       /)

nCntl = dimsizes(cntlRunNames)

year1=1991
;year2=1994
year2=1993

date1=year1*10000+0701
date2=year2*10000+0101

levels = (/ 30., 50., 70., 100., 1013.25, 1013.25 /)
averaging = (/ "tropical", "tropical", "tropical" , "tropical", "tropical", "global"/)
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
temp&time=time;-15      
temp@units = "K"
temp(:,:,:)=0.0e0
;print(temp&time)

obsfile="/glade/p/cesm/wawg/mmills/inputdata/rich_gridded_2009.nc"
fileid=addfile(obsfile,"r")
t1=floattoint(12*(year1-1958))+5
t2=floattoint(12*(year2-1958)-1)
;t2=t2+6
anomalies=fileid->anomalies ; ( time, pressure, lat, lon )
anomalies@_FillValue = -1e+30
latobs = anomalies&lat
gwobs = latRegWgt(latobs,"float",0) 
gwobs&lat=latobs
tobsTropicsAll=wgt_areaave_Wrap(anomalies(:,:,{-25:25},:),gwobs({-25:25}),1.0,0)
;print("   Max obs:"+max(tobsTropics))
;print("   Min obs:"+min(tobsTropics))
; Normalize against June 1990 - May 1991
t3=floattoint(12*(1990-1958))+5
t4=t3+11
tobsNorm=dim_avg_n_Wrap(tobsTropicsAll(t3:t4,:),0)
;new((/12,12/),"float")
;printVarSummary(tobsTropicsAll)
;print("t3="+t3+", t4="+t4)
;do mon=0,11 
;  do lev=0,11
;    tobsNorm(mon,lev)=avg(tobsTropicsAll(t3+mon:t4+mon:12,lev))
;  end do
;end do
tobsTropics=tobsTropicsAll(t1:t2,:)
;do t = t1,t2
;  tobsTropics(t-t1,:)=tobsTropicsAll(t,:)-tobsNorm(:)
;end do
;do yr = 0,year2-1-year1
;  tobsTropics(yr*12:yr*12+11,:) = tobsTropics(yr*12:yr*12+11,:) - tobsNorm(:,:)
;end do
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
gwobs&latitude=latobs
TSglob = wgt_areaave_Wrap(var_in,gwobs,1.0,0)
TStrop = wgt_areaave_Wrap(var_in(:,{-25:25},:),gwobs({-25:25}),1.0,0)
;print("TStrop:"+TStrop)

t1=(year1-1850)*12+5
t2=(year2-1850)*12-1
;t2=t2+6
TSobsGlobal = TSglob(t1:t2)
TSobsTropic = TStrop(t1:t2)
; Normalize against June 1990-May 1991
t3=(1990-1850)*12+5
t4=t3+11

TSnormGlob=avg(TSglob(t3:t4))
TSnormTrop=avg(TStrop(t3:t4))
print("TSnormTrop:"+TSnormTrop)

TSobsGlobal(:) = TSobsGlobal(:) - TSnormGlob
TSobsTropic(:) = TSobsTropic(:) - TSnormTrop
print("TSobsTropic:"+TSobsTropic)

;t4=(cntlyear2-1850)*12
;TSnormGlob=new(12,"float")
;TSnormTrop=new(12,"float")
;do mon=0,11 
;  TSnormGlob(mon)=avg(TSglob(t3+mon:t4+mon:12))
;  TSnormTrop(mon)=avg(TStrop(t3+mon:t4+mon:12))
;end do
;do yr = 0,year2-1-year1
;  TSobsGlobal(yr*12:yr*12+11) = TSobsGlobal(yr*12:yr*12+11) - TSnormGlob(0:11)
;  TSobsTropic(yr*12:yr*12+11) = TSobsTropic(yr*12:yr*12+11) - TSnormTrop(0:11)
;end do
delete(TSglob)
delete(TStrop)

datestr = systemfunc("date +%y%m%d")
datestr = datestr+"_"+systemfunc("date +%k%M")

pdffile="Tanomaly_Pinatubo_"+tag+"_c"+datestr
pdfpath="$IMAGE_PATH/"+runNames(0)+"/"
mkdir=systemfunc("mkdir -p "+pdfpath)

print("creating "+pdfpath+pdffile+".pdf")
wks = gsn_open_wks("pdf", pdfpath+pdffile)

do l = 0, nLevs - 1
  print("")
  print("Level :  " + levels(l))
  do i = 1, nRuns
    print("")
    print("  Expt Run :  " + runNames(i-1))
    print("  Cntl Run :  " + cntlRunNames(i-1))
    
    if (averaging(l).eq."global") then
      filename=runNames(i-1)+".cam.h0.globavg.nc"
      filepath=rootpath+runNames(i-1)+"/atm/proc/h0gm/"
      cntlFilename=cntlRunNames(i-1)+".cam.h0.globavg.nc"
      cntlFilepath=rootpath+cntlRunNames(i-1)+"/atm/proc/h0gm/"
    else
      filename=runNames(i-1)+".cam.h0.tropicsavg.nc"
      filepath=rootpath+runNames(i-1)+"/atm/proc/h0tm/"
      cntlFilename=cntlRunNames(i-1)+".cam.h0.tropicsavg.nc"
      cntlFilepath=rootpath+cntlRunNames(i-1)+"/atm/proc/h0tm/"
    end if
    
    fileid=addfile(filepath+filename,"r")
    cntlFileid=addfile(cntlFilepath+cntlFilename,"r")
   
    if (levels(l).eq.1013.25) then
      T      = fileid->TS
      Tcntl = cntlFileid->TS      
    else
      T      = fileid->T(:,{levels(l)})
      Tcntl = cntlFileid->T(:,{levels(l)})
    end if
    
    date=fileid->date
    d1=ind(date.eq.date1)
    d2=ind(date.eq.date2)
    
    cntldate=cntlFileid->date
    cd1=ind(cntldate.eq.date1)
    cd2=ind(cntldate.eq.date2)
    
    print("d1="+d1+" d2="+d2+" cd1="+cd1+" cd2="+cd2)
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
    if (averaging(l).eq."global") then
      temp(l, 0, :) = (/TSobsGlobal(:)/)
    else
      temp(l, 0, :) = (/TSobsTropic(:)/)
    end if
  else
;    printVarSummary(temp)
;    printVarSummary(tobsTropics)
    temp(l, 0, :) = (/tobsTropics(:,{levels(l)})/)
  end if
;  print("   Max obs:"+max(temp(l, nRuns, :)))
;  print("   Min obs:"+min(temp(l, nRuns, :)))
  res = True
  res@gsnMaximize           = True
  res@vpHeightF             = 0.5              ; change aspect ratio
  res@vpWidthF              = 1.0
  res@xyLineColors    = (/"black","blue","green4","red","orange","purple"/)
  res@xyMonoLineColor = False
  res@xyMonoDashPattern = True
  res@xyMonoLineThickness = False
  res@xyLineThicknesses = (/8.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0/)
  res@trXMinF = int2flt(year1)+5./12.
  res@trXMaxF = int2flt(year2)-31./365.
  res@tmXBMinorOn = False
;  res@tmXBMinorPerMajor = 11
  nticks=(year2-year1)*12-5
  res@tmXBMaxTicks = nticks
  res@tmXBMode= "Explicit"
  ticks=int2flt(year1)+temp&time/365.
  res@tmXBValues =ticks
  res@tmXBLabelFontHeightF = 0.015
;  XBLabels = new(nticks,"string")
;  res@tmXBLabels = (/"91","F","M","A","M","J","J","A","S","O","N","D", \
;                     "92","F","M","A","M","J","J","A","S","O","N","D", \
;                     "93","F","M","A","M","J","J","A","S","O","N","D", \
;                     "94" /)
  res@tmXBLabels = (/"J","J","A","S","O","N","D", \
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
    res@gsnCenterString = "Surface, "+averaging(l)+ \
     " average anomaly relative to June 1990 - May 1991 average"
  else
    res@gsnCenterString = levels(l)+" hPa, "+averaging(l)+\
     " average anomaly"
  end if
  
  res@gsnRightString = " "
  res@tiXAxisString = "Year"
  res@tiYAxisString = "Temperature anomaly (K)"
;  print(temp&time)
;  print(ticks)
  plot = gsn_csm_xy(wks, ticks, temp(l, :, :), res)
  
end do
print("file: "+pdfpath+pdffile+".pdf")

delete(wks)


