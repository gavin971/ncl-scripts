
histFilePath="$MODEL_PATH/waccm/run/carma-soot/bccntl02/h0/bccntl02.cam2.h0.2025-01.nc"
saveFilePath="$MODEL_PATH/waccm/analysis/nucBCemis.150Tg.4x5.wcm.nc"

histFile = addfile(histFilePath, "r")

system("/bin/rm " + saveFilePath)
saveFile=addfile(saveFilePath,"c")
;Create UNLIMITED Record Dimension
;  An UNLIMITED record dimension in a netCDF file is critical if the user 
;  plans to ever use the netCDF operators (NCO) to concatenate the file 
;  along the time/record dimension.  
filedimdef(saveFile,"time",-1,True)
;
; For a nicer looking netCDF, create a "new line" character.
; This is not necessary.
;
  nl = integertochar(10)  ; newline character
; globalAtt can be of any type. Here logical is used by convention.
;
  globalAtt             = True

  globalAtt@creation_date     = nl+\
      systemfunc("date")
      
  globalAtt@creator = nl+\
      "Mike Mills, mmills@ucar.edu"
  globalAtt@title = nl+\
      "Black carbon emission rates following US-Russia nuclear war on May 15, 2025"+nl+\
      "Integrated total BC produced over 7 days = 150 Tg"
  globalAtt@history = nl+\
       "Created by program spectrum_daily.ncl"
  globalAtt@notes = nl+\
       "Data derived from scenario described in Robock et al. Nuclear winter revisited"+nl+\
       "with a modern climate model and current nuclear arsenals: Still catastrophic"+nl+\
       "consequences. J. Geophys. Res (2007) vol. 112 pp. 13107"+nl+\
       "This scenario distributes equal mass over the 4 model layers we had between 300-150 hPa,"+nl+\
       "distributed with equal mass per gridbox over the US and Russia. Approximate rectangles"+nl+\ 
       "in the 4x5-degree lat/lon grid cover the US (75-120W, 30-54N; 10 lon x 7 lat) and"+nl+\ 
       "Russia (35-160E, 50-70N; 26 long x 6 lat). BC emissions decrease linearly to 0 by day 7:"+nl+\ 
       "  Day 1 = 42.7 Tg"+nl+\ 
       "  Day 2 = 35.8 Tg"+nl+\
       "  Day 3 = 28.6 Tg"+nl+\
       "  Day 4 = 21.5 Tg"+nl+\
       "  Day 5 = 14.3 Tg"+nl+\
       "  Day 6 = 7.1 Tg"+nl+\
       "The amounts were evenly distributed over each time step over the course of the day."+nl+\
       "These amounts are divided by the 904 gridboxes affected and 3600 s/day to derive emission"+nl+\
       "rates in kg/grid/s."

  fileattdef( saveFile, globalAtt )

lat=histFile->lat
lon=histFile->lon
lev=histFile->lev
nlat=dimsizes(lat)
nlon=dimsizes(lon)
nlev=dimsizes(lev)

;d1=2025.0*365.0+133.0
ntime=10
time=new(ntime,double)
 time@long_name="time"
 time@units = "days since 2025-05-01 00:00:00" ;
 time@calendar = "noleap"
 time(:)=(/0.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,25.0*365.0/)
 time!0="time"
print("time:");+d1)
print(time)

date=new(ntime,long)
 date(:)  =(/20250501,20250514,20250515,20250516,20250517,20250518,20250519,20250520,20250521,20500501/)
 date@long_name="date (YYYYMMDD)"
 date!0="time"
 date&time=time
print("==========================================")
print("date:")
print(date)

bctot=150 ; Tg BC emitted 
bcemisfrac=(/0.0,0.0,0.2847,0.2387,0.1907,0.1433,0.0953,0.0473,0.0,0.0/)
 bcemisfrac!0="time"
 bcemisfrac&time=time

bcemis=bctot*bcemisfrac
print("==========================================")
print("Daily BC total (Tg)")
print(bcemis)
bcemis=bcemis*1.e9/904./3600.
 bcemis!0="time"
 bcemis&time=time
print("==========================================")
print("Daily BC (kg/grid/s)")
print(bcemis)

SOOT=new((/ntime,nlev,nlat,nlon/),double)
SOOT@units="kg/grid/s"
SOOT!0="time"
SOOT!1="lev"
SOOT!2="lat"
SOOT!3="lon"
SOOT&time=time
SOOT&lev=lev
SOOT&lat=lat
SOOT&lon=lon
SOOT(:,:,:,:)=0.0

do i = 0, ntime-1
  SOOT(i,{300:150},{30:54},{240:285})=bcemis(i) ; USA
  SOOT(i,{300:150},{50:70}, {35:160})=bcemis(i) ; Russia
end do

saveFile->time=time
saveFile->date=date

saveFile->lat=histFile->lat
saveFile->lon=histFile->lon
saveFile->lev=histFile->lev
saveFile->hyai=histFile->hyai
saveFile->hybi=histFile->hybi
saveFile->hyam=histFile->hyam
saveFile->hybm=histFile->hybm
saveFile->gw=histFile->gw
saveFile->SOOT=SOOT

