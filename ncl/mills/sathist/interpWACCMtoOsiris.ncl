load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;waccmRun=getenv("caseid")

;pathList=(/"/Volumes/MillsPassport/cesm/case/",\
;           "/Volumes/Data/Models/cesm/case/",\
;           "/glade/proj2/acd/waccm/mmills/runs/",\
;           "/glade/proj3/CESM/wawg/runs/"/)

;i=0
;filePresent=False           
;do while(i.lt.dimsizes(pathList).and..not.(filePresent))
;  rootdir=pathList(i)+waccmRun+"/"
;  filePresent=isfilepresent(rootdir)
;  print(i+" "+rootdir+" "+filePresent)
;  i=i+1
;end do

;if (.not.(filePresent)) then
;  print("file not found")
;  exit
;end if

hspath=getenv("hspath")
hsfile=getenv("hsfile")
outfile=getenv("outfile")

;hsfile=rootdir+"hs/"+waccmRun+".cam.hs.02.nc"
hsfile=hspath+"/"+hsfile
print("adding "+hsfile)
hs=addfile(hsfile,"r")

osirisfile="$DATA_PATH/Osiris/osiris_hs_nabro_20110601_20110822.nc"
print("adding "+osirisfile)
osiris=addfile(osirisfile,"r")

lev_m=osiris->lev ; Osiris levels (lev=45) in meters. We will interpolate WACCM output to these levels

isOsiris=ind(hs->instr_num.eq.1)
print("interpolating "+dimsizes(isOsiris)+" Osiris columns")

z3_waccm=hs->Z3(isOsiris,:) ; WACCM geopotential height (ncol=16402, lev=88) in meters.

ext_waccm=hs->PURSULEX750(isOsiris,:) ; WACCM extinction (ncol=16402, lev=88) in /km. 

; Flip WACCM levels and interpolate to Osiris levels.
ext_waccm_osiris_levs = linint1_Wrap (z3_waccm(:,::-1), ext_waccm(:,::-1), False, lev_m, 0)

;outfilename=rootdir+"hs/"+waccmRun+".cam.hs.OsirisLevs.nc"
outfilename=hspath+"/"+outfile
print("creating "+outfilename)
fileout=addfile(outfilename,"c")
print("adding PURSULEX750: "+dimsizes(ext_waccm_osiris_levs))
fileout->PURSULEX750=ext_waccm_osiris_levs
print("adding date")
fileout->date=hs->date(isOsiris)
print("adding datesec")
fileout->datesec=hs->datesec(isOsiris)
print("adding instr_lat")
fileout->instr_lat=hs->instr_lat(isOsiris)
print("adding instr_lon")
fileout->instr_lon=hs->instr_lon(isOsiris)
print("adding obs_date")
fileout->obs_date=hs->obs_date(isOsiris)
print("adding obs_time")
fileout->obs_time=hs->obs_time(isOsiris)
print("adding orbit_num")
fileout->orbit_num=hs->orbit_num(isOsiris)
print("adding prof_num")
fileout->prof_num=hs->prof_num(isOsiris)
print("adding instr_num")
fileout->instr_num=hs->instr_num(isOsiris)
print("adding julian")
fileout->julian=hs->julian(isOsiris)
print("adding doy")
fileout->doy=hs->doy(isOsiris)
print("script complete")
