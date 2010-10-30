load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("calcZdel")
function calcZdel(a:file)

local z3,zdel
begin
 
  print("Reading geopotential heights...")
  z3 = a->Z3
  nlev=dimsizes(z3&lev)
  
  print("Calculating zdel...")
  zdel = z3
  zdel(:,0:nlev-2,:,:)=z3(:,0:nlev-2,:,:)-z3(:,1:nlev-1,:,:)
    
  delete(z3)
  delete(nlev)
  
  return(zdel)
    
end

undef("calcCol_conc")
function calcCol_conc(zdel:float, rho:float)

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

  levcol!0 = "time"
  levcol!1 = "lev"
  levcol!2 = "lat"
  levcol!3 = "lon"

  print("Calculating column ...")
;  column = dim_sum_Wrap(levcol(time|:, lat|:, lev|:))
  column = dim_sum_Wrap(levcol(time|:, lat|:, lon|:, lev|:))
  delete(levcol)
  column!0 = "time"
  column&time = rho&time
  column!1 = "lat"
  column&lat = rho&lat
  column!2 = "lon"
  column&lon = rho&lon
  print("dimsizes in calcColZavg_conc:"+dimsizes(column))
  
  return(column)
end

undef("saveSootCol")
procedure saveSootCol(historyFile:string, rootPath:string)

local historyPattern,historyNames,history,T,Tmean,pi,pm,press,temp,tmean

begin

  pi=acos(0.0)*2.0
  
  a = addfile(rootPath+historyFile,"w")

  print("  reading soot...")
  SOOTMD=a->SOOTMD  ; mass density (g/cm3)
  nlat=dimsizes(SOOTMD&lat)
  nlon=dimsizes(SOOTMD&lon)

  print("  calculating zdel...")
  zdel=calcZdel(a)  ; delta z (m)
  zdel=zdel*100.        ; m -> cm

  print("  calculating soot column...")
  SOOTMD_col=new((/1,nlat,nlon/),float)
  SOOTMD_col!0 = "time"
  SOOTMD_col&time = SOOTMD&time
  SOOTMD_col!1 = "lat"
  SOOTMD_col&lat = SOOTMD&lat
  SOOTMD_col!2 = "lon"
  SOOTMD_col&lon = SOOTMD&lon
  
  SOOTMD_col   =calcCol_conc(zdel, SOOTMD) ; g/cm2
  SOOTMD_col   =SOOTMD_col*1e13 ; Convert to ng/m2
  SOOTMD_col@units="ng/m2"

  print("  saving soot column...")
  a->SOOTMD_col=SOOTMD_col
  delete(SOOTMD_col)
  
end

;historyFile="vgsoot08_ANN_climo.nc"
;historyFile="vgsoot08_DJF_climo.nc"
historyFile="vgsoot08_JJA_climo.nc"
rootPath="/Volumes/Data/Models/waccm/run/carma-soot/vgsoot08/climo/2035-2064/"

saveSootCol(historyFile, rootPath)

