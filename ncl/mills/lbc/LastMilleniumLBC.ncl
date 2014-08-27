; Extend WACCM LBC file back to 850AD

rootpath="/glade/p/work/mmills/inputdata/lb/"
LM_ann_Filename="ghg_pmip3_850-2007_annual_c100517.v2.nc"
; vars: CO2, N2O, CH4, f11, f12
; date: 8500701 to 20070701 annual

WACCM_monthly_Filename="LBC_1765-2005_1.9x2.5_CMIP5_za_c111110.nc"
; vars: CFCL3_LBC, CF2CL2_LBC, CFC113_LBC, CCL4_LBC, CH3CCL3_LBC,
; HCFC22_LBC, CH3CL_LBC, CH3BR_LBC, CF2CLBR_LBC, CF3BR_LBC, CH4_LBC,
; CH4_LBC_mean, N2O_LBC, CO2_LBC, CO2_LBC_mean, H2_LBC, TOT_CLOY,
; TOT_BROY, OCS_LBC
; date: 17650101 to 20060101, monthly

LM_ann_File=addfile(rootpath+LM_ann_Filename,"r")
ind1765=ind(LM_ann_File->date.eq.17650701)
print("ind1765:"+ind1765)

WACCM_monthly_File=addfile(rootpath+WACCM_monthly_Filename,"r")
timeY1=WACCM_monthly_File->time(0:11)
print(timeY1)
dateY1=WACCM_monthly_File->date(0:11)
dateY1=dateY1-17650000
print(dateY1)

creation_date = systemfunc("date +%y%m%d")
outFilename="LBC_850-1764_1.9x2.5_CMIP5_za_c"+creation_date+".nc"
rm=systemfunc("rm "+rootpath+outFilename)
outFile=addfile(rootpath+outFilename,"c")
;
; For a nicer looking netCDF, create a "new line" character.
; This is not necessary.
;
  nl = integertochar(10)  ; newline character
;
; Define global attributes.
;
; globalAtt can be of any type. Here logical is used by convention.
;
  globalAtt             = True

  globalAtt@creation_date     = nl+\
      systemfunc("date")
      
  globalAtt@creator = nl+\
      "Mike Mills, mmills@ucar.edu"
  globalAtt@Title = nl+\
      "Time dependent fixed vmr surface boundary conditions"
  globalAtt@history = nl+\
      "created with program LastMilleiumLBC.ncl"+nl+\
      "extended file "+WACCM_monthly_Filename+nl+\
      "using data from "+LM_ann_Filename
      
  fileattdef( outFile, globalAtt )

lat=WACCM_monthly_File->lat
lon=WACCM_monthly_File->lon
nYears=1765-850
nMonths=nYears*12
nLats=dimsizes(lat)

print("creating time...")
time=new(nMonths,"float")
time@long_name="time"
time@units="days since 850-01-01 00:00:00"
time@calendar="noleap"
time@bounds="time_bnds"
time!0="time"

print("creating date...")
date=new(nMonths,"integer")
date@long_name="current date (YYYYMMDD)"
date!0="time"
date&time=time

do y=0,nYears-1
  year=850+y
  offset=365.d0*y
  time(y*12:y*12+11)=doubletofloat(offset)+timeY1
  date(y*12:y*12+11)=year*10000+dateY1
;  print(y+" "+year+" "+time(y*12:y*12+11)+" "+date(y*12:y*12+11))
end do
print(time(0:23))
; Create an UNLIMITED record dimension in the output netCDF file.  This is critical if 
;  the user plans to ever use NCO to concatenate the file along the time/record dimension.
filedimdef(outFile,"time",-1,True)
 
outFile->lat=lat
;outFile->time=time
outFile->date=date

vars2Da = (/ "CFCL3_LBC", "CF2CL2_LBC", "CH4_LBC", "N2O_LBC", "CO2_LBC" /)
vars_a = (/ "f11", "f12", "CH4", "N2O", "CO2" /)
vars2Db = (/ "CFC113_LBC", "CCL4_LBC", "CH3CCL3_LBC", "HCFC22_LBC", "CH3CL_LBC", \
  "CH3BR_LBC", "CF2CLBR_LBC", "CF3BR_LBC", "H2_LBC", "TOT_CLOY", "TOT_BROY", "OCS_LBC" /)
vars1D = (/ "CH4_LBC_mean", "CO2_LBC_mean" /)

print("creating variables...")
do i=0,dimsizes(vars2Da)-1
  print("...creating "+vars2Da(i))
  filevardef(outFile, vars2Da(i), typeof(WACCM_monthly_File->$vars2Da(i)$), (/"time","lat"/))
  filevarattdef(outFile, vars2Da(i), WACCM_monthly_File->$vars2Da(i)$)
  tmp=new((/nMonths,nLats/),"float")
  tmp!0="time"
  tmp!1="lat"
  tmp2=LM_ann_File->$vars_a$
  tmp3=WACCM_monthly_File->$vars2Da(i)$(0:11,:)
  do y=0,nYears-1
    year=850+y
    scalefactor=tmp2(y)/tmp2(ind1765)
    tmp(y*12:y*12+11,:)=scalefactor*tmp3
  end do
  outFile->$vars2Da(i)$=(/tmp/)
  delete(tmp)
  delete(tmp2)
  delete(tmp3)
end do

do i=0,dimsizes(vars2Db)-1
  print("...creating "+vars2Db(i))
  filevardef(outFile, vars2Db(i), typeof(WACCM_monthly_File->$vars2Db(i)$), (/"time","lat"/))
  filevarattdef(outFile, vars2Db(i), WACCM_monthly_File->$vars2Db(i)$)
  tmp=new((/nMonths,nLats/),"float")
  tmp!0="time"
  tmp!1="lat"
  tmp3=WACCM_monthly_File->$vars2Db(i)$(0:11,:)
  do y=0,nYears-1
    tmp(y*12:y*12+11,:)=tmp3
  end do
  outFile->$vars2Db(i)$=(/tmp/)
  delete(tmp)
  delete(tmp3)
end do

do i=0,dimsizes(vars1D)-1
  print("...creating "+vars1D(i))
  filevardef(outFile, vars1D(i), typeof(WACCM_monthly_File->$vars1D(i)$), (/"time"/))
  filevarattdef(outFile, vars1D(i), WACCM_monthly_File->$vars1D(i)$)
  tmp=new((/nMonths/),"float")
  tmp!0="time"
  tmp2=LM_ann_File->$vars_a$
  tmp3=WACCM_monthly_File->$vars1D(i)$(0:11)
  do y=0,nYears-1
    year=850+y
    scalefactor=tmp2(y)/tmp2(ind1765)
    tmp(y*12:y*12+11)=scalefactor*tmp3
  end do
  outFile->$vars1D(i)$=(/tmp/)
  delete(tmp)
  delete(tmp2)
  delete(tmp3)
end do

outFile->time=(/time/)

cmd="cp "+rootpath+WACCM_monthly_Filename+" "+rootpath+"tmp.nc"
print(cmd)
cp=systemfunc(cmd)
print("addfile tmp.nc")
ncf_id=addfile(rootpath+"tmp.nc","w")
adjfac=(1765.d0-850.d0)*365.d0
ncf_id->time=ncf_id->time+doubletofloat(adjfac)
ncf_id->time@units="days since 850-01-01 00:00:00"

cmd="ncrcat "+rootpath+outFilename+" "+rootpath+"tmp.nc  "+rootpath+"LBC_850-2005_1.9x2.5_CMIP5_za_c"+creation_date+".nc"
print(cmd)
ncrcat=systemfunc(cmd)

cmd="rm "+rootpath+"tmp.nc"
print(cmd)
rm=systemfunc(cmd)

