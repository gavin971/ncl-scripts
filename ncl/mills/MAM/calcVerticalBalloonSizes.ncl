load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

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

pi      = 4.0 * atan(1.0) 
sqrt2pi=sqrt(2*pi)

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

col1=strlen(caseid+".cam.h0.")
col2=col1+3
year=stringtointeger(str_get_cols(historyNames,col1,col2))
yi=ind(year.ge.startyear.and.year.le.endyear)
historyNames2=historyNames(yi)
print("historyNames2:"+dimsizes(historyNames2))

history = addfiles(rootpath+caseid+"/atm/proc/h0mam/"+historyNames2, "r")
ListSetType(history, "cat")

if isfilevar(history(0),"sigmag") then
  print("reading sigmag")
  sigmag=history[0]->sigmag
  print(sigmag)
else
  sigmag=(/1.6, 1.6, 1.2/) ; geometric standard deviations for 3 MAM3 modes
  print("sigmag not found. Using standard mode widths:")
  print(sigmag)
end if

print("reading date")
date = history[:]->date
year = floattointeger(yyyymmdd_to_yyyyfrac(date,0.0)-1./24.)
;di=ind(year.ge.startyear.and.year.le.endyear)

print("reading dgnumwet1") ; (time, lev, lat, lon)
dgnumwet1 = history[:]->dgnumwet1(:,:,{41},{105})
;printVarSummary(dgnumwet1)
print("reading dgnumwet2")
dgnumwet2 = history[:]->dgnumwet2(:,:,{41},{105})
;printVarSummary(dgnumwet2)
print("reading dgnumwet3")
dgnumwet3 = history[:]->dgnumwet3(:,:,{41},{105})
;printVarSummary(dgnumwet3)
print("reading num_a1")
num_a1 = history[:]->num_a1(:,:,{41},{105})
;printVarSummary(num_a1)
print("reading num_a2")
num_a2 = history[:]->num_a2(:,:,{41},{105})
;printVarSummary(num_a2)
print("reading num_a3")
num_a3 = history[:]->num_a3(:,:,{41},{105})
;printVarSummary(num_a3)


bin=(/0.,0.15,0.25,0.30,0.49,0.62,0.78,1.08/)
bin@units="um"
bin@long_name="minimum radius"
nbins=dimsizes(bin)

