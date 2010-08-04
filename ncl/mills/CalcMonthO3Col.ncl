load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("calcPress")
function calcPress(history:list, historyNames[*]:string)

local pi,pm,press,wgt
begin
  pi=acos(0.0)*2.0

  print("Reading surface pressure...")
  ps = addfiles_GetVar(history, historyNames, "PS")
  p0 = history[0]->P0
  hyam = history[0]->hyam
  hybm = history[0]->hybm
  time = addfiles_GetVar(history, historyNames, "time")
  lev = history[0]->lev
  lat = history[0]->lat
  lon = history[0]->lon
  
  print("Calculating pressure at midpoints...")
  pm = pres_hybrid_ccm(ps,p0,hyam,hybm)
  print(dimsizes(pm))
  pm!0 = "time"
  pm!1 = "lev"
  pm&lev = lev
  pm!2 = "lat"
  pm&lat = lat
  pm!3 = "lon"
  pm&lon = lon
  
  delete(p0)
  delete(ps)
  delete(hyam)
  delete(hybm)
  
  print("averaging pressure over time")
  press=dim_avg_Wrap(pm(lev|:, lat|:, lon|:, time|:))
  
  delete(pm)
  
  print("averaging pressure zonally")
  pm=dim_avg_Wrap(press(lev|:, lat|:, lon|:))
  
  delete(press)
  
  wgt=dble2flt(cos(lat*pi/180.))
  
  print("averaging pressure meridionally")
  pressure=dim_avg_wgt_Wrap(pm(lev|:, lat|:),wgt,0)
  
  delete(pm)
  
  pressure=pressure/100.
  
  pressure@units = "hPa"

  return(pressure)
end

undef("calcPdel")
function calcPdel(history:list, historyNames[*]:string)

local g,Mair,partPress,column,vmr,pdel
begin
 
  print("Reading surface pressure...")
  ps = addfiles_GetVar(history, historyNames, "PS")
  p0 = history[0]->P0
  hyai = history[0]->hyai
  hybi = history[0]->hybi
  time = addfiles_GetVar(history, historyNames, "time")
  lev = history[0]->lev
  lat = history[0]->lat
  lon = history[0]->lon
  
  print("Calculating pressure at interfaces...")
  pm = pres_hybrid_ccm(ps,p0,hyai,hybi)
    
  delete(p0)
  delete(ps)
  
  delete(hyai)
  delete(hybi)
  
  print("Calculating pressure change...")
  pdel=pm(:,1:dimsizes(lev),:,:)-pm(:,0:dimsizes(lev)-1,:,:)
  
  delete(pm)
  
  return(pdel)
    
end

undef("calcZ3")
function calcZ3(history:list, historyNames[*]:string)

local z3,zdel
begin
 
  print("Reading geopotential heights...")
  z3 = addfiles_GetVar(history, historyNames, "Z3")
  nlev=dimsizes(z3&lev)
  
  z0 = z3(:,nlev-1,:,:)

  delete(z3)
  delete(nlev)
  
  return(z0)
    
end

undef("calcCol_vmr")
function calcCol_vmr(history:list, historyNames[*]:string, var:string, pdel:float)

local g,Mair,partPress,column,vmr,columnAv
begin
  ; Some constants we will need:
  g       = 980.6               ; cm s^-2 
  Mair    = 28.964/6.022e23     ; g/molec
  pi=acos(0.0)*2.0
  print("calculating column")

  ; volume mixing ratios:
  print("Reading vmr...")
  vmr=addfiles_GetVar(history, historyNames, var)
  lat=vmr&lat
  lon=vmr&lon
;  vmr=vmr/divfac

  ; dP=rho*g*dh=[M]*Mair*g*dh
  ; [M]*dh=dP/(Mair*g)
  partPress=pdel*vmr*10./(Mair*g)  ; pdel(Pa)*10(dynes/cm2/Pa)
  partPress!0 = "time"
  partPress!1 = "lev"
  partPress&lev = vmr&lev
  partPress!2 = "lat"
  partPress&lat = lat
  partPress!3 = "lon"
  partPress&lon = vmr&lon
  
  delete(vmr)
  
  print("Calculating column ...")
  column = dim_sum_Wrap(partPress(time|:, lat|:, lon|:, lev|:))
  column = column/2.687e16 ; Convert to Dobson Units
  column!0 = "time"
  column!1 = "lat"
  column&lat = lat
  column!2 = "lon"
  column&lon = lon
  column@units = "Dobson Units"
  delete(lat)

  return(column)
end

undef("saveMeanFiles")
procedure saveMeanFiles(historyRun:string, historyPath:string, subDir:string)

local historyPattern1,historyNames1,history1,historyPattern2,historyNames2, \
  history2,T,Tmean,pi,pm,press,temp,tmean,saveFile1,saveFile2

begin
  pattern="h0"
  historyPath2=historyPath+historyRun+"/"
  
  saveFile = addfile(historyPath2+historyRun+".O3Col.h0.nc", "c")
  historyPattern=historyPath2+subDir+"/*.cam2."+pattern+"*.nc"
  
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  
  ; Ozone:
  print("reading ozone...")
  press=calcPress(history, historyNames)
  pdel=calcPdel(history, historyNames)
  O3_col=calcCol_vmr(history, historyNames, "O3", pdel)
  Z0=calcZ3(history, historyNames)
  
  print("saving ozone...")
  saveFile->O3_col=O3_col
  print("saving Z0...")
  saveFile->Z0=Z0
  print("Done!")

end


historyPath="/Volumes/Data/Models/waccm/run/carma-soot/"
historyRun="Soot5Tg"
subDir="h0"

;saveMeanFiles(historyRun,historyPath,subDir)

historyRun="PassiveSoot"

saveMeanFiles(historyRun,historyPath,subDir)

