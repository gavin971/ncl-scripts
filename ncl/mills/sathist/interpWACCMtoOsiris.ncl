load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

waccmRun=getenv("caseid")

pathList=(/"/Volumes/MillsPassport/cesm/case/",\
           "/Volumes/Data/Models/cesm/case/"/)

i=0
filePresent=False           
do while(i.lt.dimsizes(pathList).and..not.(filePresent))
  rootdir=pathList(i)+waccmRun+"/"
  filePresent=isfilepresent(rootdir)
  print(i+" "+rootdir+" "+filePresent)
  i=i+1
end do

if (.not.(filePresent)) then
  print("file not found")
  exit
end if

hs=addfile(rootdir+"hs/"+waccmRun+".cam.hs.02.nc","r")

osiris=addfile("$DATA_PATH/Osiris/osiris_hs_nabro_20110601_20110822.nc","r")

lev_m=osiris->lev ; Osiris levels (lev=45) in meters. We will interpolate WACCM output to these levels

isOsiris=ind(hs->instr_num.eq.1)
print(dimsizes(isOsiris))

z3_waccm=hs->Z3(isOsiris,:) ; WACCM geopotential height (ncol=16402, lev=88) in meters.

ext_waccm=hs->PURSULEX750(isOsiris,:) ; WACCM extinction (ncol=16402, lev=88) in /km. 

; Flip WACCM levels and interpolate to Osiris levels.
ext_waccm_osiris_levs = linint1_Wrap (z3_waccm(:,::-1), ext_waccm(:,::-1), False, lev_m, 0)

outfile=addfile(rootdir+"hs/"+waccmRun+".cam.hs.OsirisLevs.nc","c")
outfile->PURSULEX750=ext_waccm_osiris_levs
outfile->date=hs->date(isOsiris)
outfile->datesec=hs->datesec(isOsiris)
outfile->instr_lat=hs->instr_lat(isOsiris)
outfile->instr_lon=hs->instr_lon(isOsiris)
outfile->obs_date=hs->obs_date(isOsiris)
outfile->obs_time=hs->obs_time(isOsiris)
outfile->orbit_num=hs->orbit_num(isOsiris)
outfile->prof_num=hs->prof_num(isOsiris)
outfile->instr_num=hs->instr_num(isOsiris)
outfile->julian=hs->julian(isOsiris)
outfile->doy=hs->doy(isOsiris)
