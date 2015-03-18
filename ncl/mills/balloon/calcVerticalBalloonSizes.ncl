load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

printflag=.not.ismissing(getenv("printflag"))

caseid=getenv("caseid")

if (ismissing(caseid)) then
  print("FATAL: caseid environment variable not set")
  exit
else
  print("caseid: "+caseid)
end if

startyear=stringtoint(getenv("startyear"))
endyear=stringtoint(getenv("endyear"))

if (ismissing(startyear)) then
  startyear=1999
end if  

if (ismissing(endyear)) then
  endyear=2008
end if

if (startyear.gt.endyear) then
  print("FATAL: startyear = "+startyear+" > endyear = "+endyear)
  exit
end if

print("startyear = "+startyear+", endyear = "+endyear)

if (startyear.eq.endyear) then
  rangestr=startyear
else
  rangestr=startyear+"-"+endyear
end if

sqrt2=sqrt(2.0)

rootpath = "/glade/scratch/mmills/archive/"

exists =isfilepresent(rootpath+caseid+"/atm/proc/h0mam/")
if (.not.exists) then
  print("Missing: "+rootpath+caseid+"/atm/proc/h0mam/")
  rootpath="/glade/p/cesm/wawg_dev/mmills/archive/"
end if

exists =isfilepresent(rootpath+caseid+"/atm/proc/h0mam/")
if (.not.exists) then
  print("Missing: "+rootpath+caseid+"/atm/proc/h0mam/")
  print("FATAL: check caseid="+caseid)
  exit
end if
print("Rootpath: "+rootpath+caseid+"/atm/proc/h0mam/")

historyNames=systemfunc("ls "+rootpath+caseid+"/atm/proc/h0mam/")
print("historyNames:"+dimsizes(historyNames))

; select only files between startyear and endyear
col1=strlen(caseid+".cam.h0.")
col2=col1+3
year=stringtointeger(str_get_cols(historyNames,col1,col2))
yi=ind(year.ge.startyear.and.year.le.endyear)
historyNames2=historyNames(yi)
print("historyNames2:"+dimsizes(historyNames2))

history = addfiles(rootpath+caseid+"/atm/proc/h0mam/"+historyNames2, "r")
ListSetType(history, "cat")

if isfilevar(history[0],"sigmag") then
  print("reading sigmag")
  sigmag=history[0]->sigmag
  print(sigmag)
else
  sigmag=(/1.6, 1.6, 1.2/) ; geometric standard deviations for 3 MAM3 modes
  print("sigmag not found. Using standard mode widths:")
  print(sigmag)
end if
logsigmag=log(sigmag)

;print("reading date")
;date = history[:]->date
;year = floattointeger(yyyymmdd_to_yyyyfrac(date,0.0)-1./24.)
;di=ind(year.ge.startyear.and.year.le.endyear) 

; Laramie (41N, 105W)
lat1=41
lon1=360-105

print("reading dgnumwet1") ; (time, lev, lat, lon)
dgnumwet1 = history[:]->dgnumwet1(:,:,{lat1},{lon1}); Aerosol mode wet diameter in meters
;printVarSummary(dgnumwet1)
print("reading dgnumwet2")
dgnumwet2 = history[:]->dgnumwet2(:,:,{lat1},{lon1})
;printVarSummary(dgnumwet2)
print("reading dgnumwet3")
dgnumwet3 = history[:]->dgnumwet3(:,:,{lat1},{lon1})
;printVarSummary(dgnumwet3)

print("reading dgnd_a01") ; (time, lev, lat, lon)
dgnd_a01 = history[:]->dgnd_a01(:,:,{lat1},{lon1}); Aerosol mode dry diameter in meters
;printVarSummary(dgnd_a01)
print("reading dgnd_a02")
dgnd_a02 = history[:]->dgnd_a02(:,:,{lat1},{lon1})
;printVarSummary(dgnd_a02)
print("reading dgnd_a03")
dgnd_a03 = history[:]->dgnd_a03(:,:,{lat1},{lon1})
;printVarSummary(dgnd_a03)

print("reading num_a1")
num_a1 = history[:]->num_a1(:,:,{lat1},{lon1}) ; aerosol number concentration in /kg air
;printVarSummary(num_a1)
print("reading num_a2")
num_a2 = history[:]->num_a2(:,:,{lat1},{lon1})
;printVarSummary(num_a2)
print("reading num_a3")
num_a3 = history[:]->num_a3(:,:,{lat1},{lon1})
;printVarSummary(num_a3)

print("reading arhomass")
arhomass = history[:]->arhomass(:,:,{lat1},{lon1}) ; air density in g/cm3
print("reading Z3")
Z3 = history[:]->Z3(:,:,{lat1},{lon1}) ; geopotential height (m)

diminq=dimsizes(dgnumwet1)
ntimes=diminq(0)
nlev=diminq(1)

; convert #/kg air to #/cm3
num_a1=num_a1*arhomass/1.e3
num_a2=num_a2*arhomass/1.e3
num_a3=num_a3*arhomass/1.e3

bin=(/0.01,0.15,0.25,0.30,0.49,0.62,0.78,1.08/)
bin@units="um"
bin@long_name="minimum radius"
nbins=dimsizes(bin)

; diameter = 2* radius, 1e-6 m/um
dbin=bin*2e-6

nd=new((/nlev,nbins,ntimes/),"float")
nd1=nd
nd2=nd
nd3=nd
print("calculating number densities")
do t =0,ntimes-1
  t1=t+1
  print("time "+t1+" of "+ntimes)
  do k =0,nlev-1
    do b=0,nbins-1
      e1 = erf(log(dbin(b)/dgnumwet1(t,k)/sqrt2/logsigmag(0)))
      e2 = erf(log(dbin(b)/dgnumwet2(t,k)/sqrt2/logsigmag(1)))
      e3 = erf(log(dbin(b)/dgnumwet3(t,k)/sqrt2/logsigmag(2)))
      nd1(k,b,t) = num_a1(t,k)*(0.5-0.5*e1)
      nd2(k,b,t) = num_a2(t,k)*(0.5-0.5*e2)
      nd3(k,b,t) = num_a3(t,k)*(0.5-0.5*e3)
      nd(k,b,t) = nd1(k,b,t) + nd2(k,b,t) + nd3(k,b,t)
      if (k.eq.45..and.printflag) then
        print("  ")
        print("lev="+dgnumwet1&lev(k)+", alt="+Z3(t,k)+" m")
        print("calculating population with radius larger than "+bin(b)+" um")
        print("diameter larger than "+dbin(b)+" m")
        print("dgnumwet (m), mode1:"+ dgnumwet1(t,k)+" mode2:"+dgnumwet2(t,k)+" mode3:"+dgnumwet3(t,k))
        print("dgn  dry (m), mode1:"+ dgnd_a01(t,k)+" mode2:"+dgnd_a02(t,k)+" mode3:"+dgnd_a03(t,k))
        print("sigmag:"+sigmag(0)+" "+sigmag(1)+" "+sigmag(2))
        print("lnsigmag:"+logsigmag(0)+" "+logsigmag(1)+" "+logsigmag(2))
        print("sqrt2:"+sqrt2)
        print("number density (/cm3), mode1:"+num_a1(t,k)+" mode2:"+num_a2(t,k)+" mode3:"+num_a3(t,k))
        print("erf, mode1:"+e1+" mode2:"+e2+" mode3:"+e3)
        print("population above threshold, mode1:"+nd1(k,b,t)+" mode2="+nd2(k,b,t)+" mode3="+nd3(k,b,t))
        print("total population="+nd(k,b,t))
      end if
    end do
  end do
end do

ndbin=dim_avg(nd)
ndbin!0="lev"
ndbin!1="bin"
ndbin&lev=dgnumwet1&lev
ndbin&bin=bin
ndbin@long_name="aerosol number density larger than bin threshold"
ndbin@units="cm-3"

nd1bin=ndbin ; get metadata
nd1bin=dim_avg(nd1)
nd1bin@long_name="aerosol mode 1 number density larger than bin threshold"

nd2bin=ndbin ; get metadata
nd2bin=dim_avg(nd2)
nd2bin@long_name="aerosol mode 2 number density larger than bin threshold"

nd3bin=ndbin ; get metadata
nd3bin=dim_avg(nd3)
nd3bin@long_name="aerosol mode 3 number density larger than bin threshold"

Z3av=dim_avg_Wrap(Z3(lev|:,time|:))

print(sprintf("%3.2e",ndbin&lev)+" "+sprintf("%3.2e",ndbin(:,0))+" "+\
      sprintf("%3.2e",ndbin(:,1))+" "+sprintf("%3.2e",ndbin(:,2))+" "+sprintf("%3.2e",ndbin(:,3)))

outfile=rootpath+caseid+"/atm/proc/"+caseid+".US_Laramie_41N_105W_"+rangestr+".nc"
print("rm "+outfile)
a=systemfunc("rm "+outfile)
print("creating "+outfile)
ncid=addfile(outfile,"c")
print("saving Z3")
ncid->Z3=Z3av
print("saving ndbin")
ncid->ndbin=ndbin
print("saving nd1bin")
ncid->nd1bin=nd1bin
print("saving nd2bin")
ncid->nd2bin=nd2bin
print("saving nd3bin")
ncid->nd3bin=nd3bin

