load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

rootdir = "/data4/mmills/f2000.e10r02.2deg.waccm.005/"

; h1: CO2,H,NO,O,O3,P0,date,hyam,hybm,lat,lev,lon,time
a=addfile(rootdir+"f2000.e10r02.2deg.waccm.005.cam2.h1zm.varscat.nc","r")

; h2: QRS_TOT,P0,date,hyam,hybm,lat,lev,lon,time
b=addfile(rootdir+"f2000.e10r02.2deg.waccm.005.cam2.h2zm.varscat.nc","r")

time=a->time
date=a->date
CO2=a->CO2
H=a->H
NO=a->NO
O=a->O
O3=a->O3
P0=a->P0
hyam=a->hyam
hybm=a->hybm
lat=a->lat
lev=a->lev
lon=a->lon

QRS_TOT=b->QRS_TOT

delete(a)
delete(b)

modtime=mod(time,365)
modtime2=mod(time-5,365)
modtime=where(modtime.eq.0, 365, modtime)
modtime2=where(modtime2.eq.0, 365, modtime2)

date_out=date(0:72)
time_out=time(0:72)
CO2_out=CO2(0:72,:,:)
H_out=H(0:72,:,:)
NO_out=NO(0:72,:,:)
O_out=O(0:72,:,:)
O3_out=O3(0:72,:,:)
QRS_TOT_out=QRS_TOT(0:72,:,:)

do i = 0,35
  j=doubletoint(modtime(i)/5-1)
  k=doubletoint(modtime2(i)/5-1)
  print(i+" "+j+" "+k+" "+modtime(i)+" "+modtime2(i))
  tind=ind(modtime.eq.modtime(i))
  date_out(j)=monthday(1,doubletoint(modtime(i)))
  time_out(j)=modtime(i)
  CO2_out(j,:,:)=dim_avg_Wrap(CO2(lev|:,lat|:,time|tind))
  H_out  (j,:,:)=dim_avg_Wrap(H  (lev|:,lat|:,time|tind))
  NO_out (j,:,:)=dim_avg_Wrap(NO (lev|:,lat|:,time|tind))
  O_out  (j,:,:)=dim_avg_Wrap(O  (lev|:,lat|:,time|tind))
  O3_out (j,:,:)=dim_avg_Wrap(O3 (lev|:,lat|:,time|tind))
  QRS_TOT_out(k,:,:)=dim_avg_Wrap(QRS_TOT(lev|:,lat|:,time|tind))
end do

do i = 37,72
  j=doubletoint(modtime(i)/5-1)
  k=doubletoint(modtime2(i)/5-1)
  print(i+" "+j+" "+k+" "+modtime(i)+" "+modtime2(i))
  tind=ind(modtime.eq.modtime(i))
  date_out(j)=monthday(1,doubletoint(modtime(i)))
  time_out(j)=modtime(i)
  CO2_out(j,:,:)=dim_avg_Wrap(CO2(lev|:,lat|:,time|tind))
  H_out  (j,:,:)=dim_avg_Wrap(H  (lev|:,lat|:,time|tind))
  NO_out (j,:,:)=dim_avg_Wrap(NO (lev|:,lat|:,time|tind))
  O_out  (j,:,:)=dim_avg_Wrap(O  (lev|:,lat|:,time|tind))
  O3_out (j,:,:)=dim_avg_Wrap(O3 (lev|:,lat|:,time|tind))
  QRS_TOT_out(k,:,:)=dim_avg_Wrap(QRS_TOT(lev|:,lat|:,time|tind))  
end do

i = 36
j=doubletoint(modtime(i)/5-1)
k=doubletoint(modtime2(i)/5-1)
print(i+" "+j+" "+k+" "+modtime(i)+" "+modtime2(i))
tind=ind(modtime.eq.modtime(i))
date_out(j)=monthday(1,doubletoint(modtime(i)))
time_out(j)=modtime(i)
CO2_out(j,:,:)=dim_avg_Wrap(CO2(lev|:,lat|:,time|tind))
H_out  (j,:,:)=dim_avg_Wrap(H  (lev|:,lat|:,time|tind))
NO_out (j,:,:)=dim_avg_Wrap(NO (lev|:,lat|:,time|tind))
O_out  (j,:,:)=dim_avg_Wrap(O  (lev|:,lat|:,time|tind))
O3_out (j,:,:)=dim_avg_Wrap(O3 (lev|:,lat|:,time|tind))
QRS_TOT_out(k,:,:)=dim_avg_Wrap(QRS_TOT(lev|:,lat|:,time|tind))

date_out&time=time_out
CO2_out&time=time_out
H_out&time=time_out
NO_out&time=time_out
O_out&time=time_out
O3_out&time=time_out
QRS_TOT_out&time=time_out

c=addfile(rootdir+"ghg_forcing_2000_c121016.nc","c")

nl = integertochar(10)  ; newline character
  
globalAtt             = True

globalAtt@creation_date     = nl+\
    systemfunc("date")

globalAtt@creator = nl+\
    "Mike Mills, mmills@ucar.edu"
globalAtt@user = nl+\
    "CESM1.0 user"
globalAtt@title = nl+\
    "Driver file for SC-WACCM (specified chemistry)" ;
globalAtt@history = nl+\
    "From years 7-16 of run f2000.e10r02.2deg.waccm.005"+nl+\
    "Zonally averaged with NCO script wazmdir" +nl+\
    "ghg_forcing_2000_c121015.nc created with NCL script waccmsc-driver.ncl"
globalAtt@notes = nl+\
    "Instantaneous output every 5 days for CO2, H, NO, O, and O3."+nl+\
    "10-day average midnight GMT values for QRS_TOT centered at 5-day intervals."
      
fileattdef( c, globalAtt )

c->date=date_out
c->P0=P0
c->hyam=hyam
c->hybm=hybm
c->lat=lat
c->lev=lev
c->lon=lon
c->CO2=CO2_out
c->H=H_out
c->NO=NO_out
c->O=O_out
c->O3=O3_out
c->QRS_TOT=QRS_TOT_out

delete(c)




