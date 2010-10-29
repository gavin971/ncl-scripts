load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("calcZdelZavg")
function calcZdelZavg(a:file)

local z3,zdel
begin
 
  print("Reading geopotential heights...")
  z3 = a->Z3
  nlev=dimsizes(z3&lev)
  
  print("Calculating zdel...")
  zdel = z3
;  zdel(:,0:nlev-2,:)=z3(:,0:nlev-2,:)-z3(:,1:nlev-1,:)
  zdel(0:nlev-2,:)=z3(0:nlev-2,:)-z3(1:nlev-1,:)
    
  delete(z3)
  delete(nlev)
  
  return(zdel)
    
end

undef("calcColZavg_vmr")
function calcColZavg_vmr(pdel:float, vmr:float)

local g,Mair,partPress,column,vmr,columnAv
begin
  ; Some constants we will need:
  g       = 980.6               ; cm s^-2 
  Mair    = 28.964/6.022e23     ; g/molec
  pi=acos(0.0)*2.0

  ; dP=rho*g*dh=[M]*Mair*g*dh
  ; [M]*dh=dP/(Mair*g)
  partPress=pdel*vmr*10./(Mair*g)  ; pdel(Pa)*10(dynes/cm2/Pa)
  partPress!0 = "time"
  partPress&time = vmr&time
  partPress!1 = "lev"
  partPress&lev = vmr&lev
  partPress!2 = "lat"
  partPress&lat = vmr&lat
    
  print("  Calculating column ...")
  column = dim_sum_Wrap(partPress(time|:, lat|:, lev|:))

  return(column)
end

undef("calcColZavg_conc")
function calcColZavg_conc(zdel:float, rho:float)

local conc,levcol,column,columnAv
begin
  ; Some constants we will need:
  g       = 980.6               ; cm s^-2 
  Mair    = 28.964/6.022e23     ; g/molec
  pi=acos(0.0)*2.0
  print("calculating column")

  print("rho:"+dimsizes(rho))
  print("zdel:"+dimsizes(zdel))
  levcol=rho*zdel

;  levcol!0 = "time"
;  levcol!1 = "lev"
;  levcol!2 = "lat"

  levcol!0 = "lev"
  levcol!1 = "lat"

  print("Calculating column ...")
;  column = dim_sum_Wrap(levcol(time|:, lat|:, lev|:))
  column = dim_sum_Wrap(levcol(lat|:, lev|:))
  delete(levcol)
;  column!0 = "time"
;  column&time = rho&time
;  column!1 = "lat"
;  column&lat = rho&lat
  column!0 = "lat"
  column&lat = rho&lat
  print("dimsizes in calcColZavg_conc:"+dimsizes(column))
  
  return(column)
end

undef("saveZavgCol")
procedure saveZavgCol(historyRun:string, rootPath:string)

local historyPattern,historyNames,history,T,Tmean,pi,pm,press,temp,tmean

begin

  pi=acos(0.0)*2.0
  
;  historyPattern = rootPath+"/"+historyRun+"/h0/*.h0.*.nc"
  historyPattern = rootPath+"/"+historyRun+"/h0zm/"+historyRun+".cam2.h0zm.[01][0-9].nc"
  historyNames = systemfunc("ls " + historyPattern)
;  zmNames = str_sub_str(historyNames,"h0","h0zm")
  zmNames = historyNames
  nFiles=dimsizes(historyNames)
   
  i = 0
  do while ((i.lt.nFiles).and.(isfilepresent(zmNames(i))))
    print("Reading from "+zmNames(i))
    a = addfile(zmNames(i),"w")
    
    ; Ozone:
;    print("  reading ozone...")
    O3=a->O3
    pdel=a->pdel
    nlat=dimsizes(O3&lat)
    
;    O3_col=new((/1,nlat/),float)
;    O3_col!0 = "time"
;    O3_col&time = O3&time
;    O3_col!1 = "lat"
;    O3_col&lat = O3&lat
    
;    O3_col=calcColZavg_vmr(pdel, O3)
;    O3_col=O3_col/2.687e16 ; Convert to Dobson Units
;    O3_col@units="Dobson Units"

;    print("  saving ozone column...")
;    a->O3_col=O3_col

    print("  reading soot...")
    SOOTMD=a->SOOTMD  ; mass density (g/cm3)
    
    print("  calculating zdel...")
    zdel=calcZdelZavg(a)  ; delta z (m)
    zdel=zdel*100.        ; m -> cm

    print("  calculating soot column...")
;    SOOT01_col=new((/1,nlat/),float)
;    SOOT01_col!0 = "time"
;    SOOT01_col&time = O3&time
;    SOOT01_col!1 = "lat"
;    SOOT01_col&lat = O3&lat
    SOOTMD_col=new(nlat,float)
    SOOTMD_col!0 = "lat"
    SOOTMD_col&lat = O3&lat
    SOOTMD_col   =calcColZavg_conc(zdel, SOOTMD) ; g/cm2
    SOOTMD_col   =SOOTMD_col*1e13 ; Convert to ng/m2
    SOOTMD_col@units="ng/m2"
    
    print("  saving soot column...")
    a->SOOTMD_col=SOOTMD_col
    delete(SOOTMD_col)

    i = i + 1
  end do
  
end

undef("saveZavgPdel")
procedure saveZavgPdel(historyRun:string, rootPath:string)
 
local historyPattern,historyNames,i,a,b,ps,p0,hyai,hybi,time,lev,lat,lon,pm,pdel,pdelZavg
 
begin
;  historyPattern = rootPath+"/"+historyRun+"/h0/*.h0.*.nc"
  historyPattern = rootPath+"/"+historyRun+"/h0zm/"+historyRun+".cam2.h0zm.[01][0-9].nc"
  print(historyPattern)
  historyNames = systemfunc("ls " + historyPattern)
;  print(historyNames)
;  zmNames = str_sub_str(historyNames,"h0","h0zm")
  zmNames = historyNames
  print(zmNames)
  nFiles=dimsizes(historyNames)
  
  i = 0
  print(i+" | "+nFiles+" | "+zmNames(i))
  print(isfilepresent(zmNames(i)))
  do while ((i.lt.nFiles).and.(isfilepresent(zmNames(i))))
    print("Reading from "+historyNames(i))
    a = addfile(historyNames(i),"r")

    print("  Reading surface pressure...")
    ps = a->PS
    p0 = a->P0
    hyai = a->hyai
    hybi = a->hybi
    time = a->time
    lev = a->lev
    lat = a->lat
    lon = a->lon
    delete(a)

    print("  Calculating pressure at interfaces...")
    pm = pres_hybrid_ccm(ps,p0,hyai,hybi)

    delete(p0)
    delete(ps)

    delete(hyai)
    delete(hybi)

    print("  Calculating pressure change...")
    pdel=pm(:,1:dimsizes(lev),:,:)-pm(:,0:dimsizes(lev)-1,:,:)
    pdel!0="time"
    pdel&time=time
    pdel!1="lev"
    pdel&lev=lev
    pdel!2="lat"
    pdel&lat=lat
    pdel!3="lon"
    pdel&lon=lon

    delete(pm)

    pdelZavg=dim_avg_Wrap(pdel(time|:, lev|:, lat|:, lon|:))

    delete(pdel)

    b = addfile(zmNames(i),"w")

    b->pdel=pdelZavg

    delete(pdelZavg)
    delete(b)
    
    i = i + 1
  end do
end

;historyRun="vgcontrol02"
historyRun="vgsoot08"
;rootPath="/data5/mmills"
rootPath="/Volumes/Data/Models/waccm/run/carma-soot"

;saveZavgPdel(historyRun, rootPath)
saveZavgCol(historyRun, rootPath)

