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
  rank=dimsizes(dimsizes(ps))
  p0 = history[0]->P0
  gw = history[0]->gw
  hyam = history[0]->hyam
  hybm = history[0]->hybm
  time = addfiles_GetVar(history, historyNames, "time")
  lev = history[0]->lev
  lat = history[0]->lat
  if (rank.eq.3) then
    lon = history[0]->lon
  end if
  
  print("Calculating pressure at midpoints...")
  pm = pres_hybrid_ccm(ps,p0,hyam,hybm)
  print(dimsizes(pm))
  if (rank.eq.3) then
    pm!0 = "time"
    pm!1 = "lev"
    pm&lev = lev
    pm!2 = "lat"
    pm&lat = lat
    pm!3 = "lon"
    pm&lon = lon
  else
    pm!1 = "time"
    pm!0 = "lev"
    pm&lev = lev
    pm!2 = "lat"
    pm&lat = lat  
  end if
  
  delete(p0)
  delete(ps)
  delete(hyam)
  delete(hybm)
  
  if (rank.eq.3) then
    print("averaging pressure zonally")
    press=dim_avg_Wrap(pm(lev|:, lat|:, time|:, lon|:))
    delete(pm)
    pm=press
    delete(press)
  end if
    
  print("averaging pressure over time")
  press=dim_avg_Wrap(pm(lev|:, lat|:, time|:))
  
  delete(pm)
  
  print("averaging pressure meridionally")
  pressure=dim_avg_wgt_Wrap(press(lev|:, lat|:),gw,0)
  
  delete(press)
  
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
  rank=dimsizes(dimsizes(ps))
  p0 = history[0]->P0
  hyai = history[0]->hyai
  hybi = history[0]->hybi
  time = addfiles_GetVar(history, historyNames, "time")
  lev = history[0]->lev
  lat = history[0]->lat
  if (rank.eq.3) then
    lon = history[0]->lon
  end if
  
  print("Calculating pressure at interfaces...")
  pm = pres_hybrid_ccm(ps,p0,hyai,hybi)
    
  delete(p0)
  delete(ps)
  
  delete(hyai)
  delete(hybi)
  
  print("Calculating pressure change...")
  if (rank.eq.3) then
    pdel=pm(:,1:dimsizes(lev),:,:)-pm(:,0:dimsizes(lev)-1,:,:)
  else
    pm!1 = "time"
    pm!0 = "lev"
    pm!2 = "lat"
    pdel=pm(time|:,lev|1:dimsizes(lev),lat|:)-pm(time|:,lev|0:dimsizes(lev)-1,lat|:)
  end if
  
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
  rank=dimsizes(dimsizes(z3))
  
  if (rank.eq.4) then
    z0 = z3(:,nlev-1,:,:)
  else
    z0 = z3(:,nlev-1,:)
  end if

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
  rank=dimsizes(dimsizes(vmr))
  time=vmr&time
  lat=vmr&lat
  if (rank.eq.4) then
    lon=vmr&lon
  end if
;  vmr=vmr/divfac

  ; dP=rho*g*dh=[M]*Mair*g*dh
  ; [M]*dh=dP/(Mair*g)
  partPress=pdel*vmr*10./(Mair*g)  ; pdel(Pa)*10(dynes/cm2/Pa)
  partPress!0 = "time"
  partPress&time = time
  partPress!1 = "lev"
  partPress&lev = vmr&lev
  partPress!2 = "lat"
  partPress&lat = lat
  if (rank.eq.4) then
    partPress!3 = "lon"
    partPress&lon = vmr&lon
  end if
  
  delete(vmr)
  
  print("Calculating column ...")
  if (rank.eq.4) then
    column = dim_sum_Wrap(partPress(time|:, lat|:, lon|:, lev|:))
  else
    column = dim_sum_Wrap(partPress(time|:, lat|:, lev|:))
  end if
  column = column/2.687e16 ; Convert to Dobson Units
  column!0 = "time"
  column&time = time
  column!1 = "lat"
  column&lat = lat
  if (rank.eq.4) then
    column!2 = "lon"
    column&lon = lon
  end if
  column@units = "Dobson Units"
  delete(lat)

  return(column)
end

undef("saveMeanFiles")
procedure saveMeanFiles(historyRun:string, historyPath:string, subDir:string)

local historyPattern1,historyNames1,history1,historyPattern2,historyNames2, \
  history2,T,Tmean,pi,pm,press,temp,tmean,saveFile1,saveFile2

begin
  pattern=subDir
  historyPath2=historyPath+historyRun+"/"
  analysisPath=historyPath2+"proc/"
  a = systemfunc("mkdir -p "+analysisPath)
  
  historyPattern=historyPath2+subDir+"/*.cam2."+pattern+"*.nc"
  
  historyNames = systemfunc("ls " + historyPattern)
  print("opening history files...")
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  
  ; Ozone:
  print("reading ozone...")
  press=calcPress(history, historyNames)
  pdel=calcPdel(history, historyNames)
  O3_col=calcCol_vmr(history, historyNames, "O3", pdel)
  O3_col@long_name="O3 column"
  Z0=calcZ3(history, historyNames)
  print("reading date...")
  date=addfiles_GetVar(history, historyNames, "date")
  
  newfile=analysisPath+historyRun+".O3Col."+pattern+".nc"
  exists = isfilepresent(newfile)
  if (exists) then
    print("Deleting existing file: "+newfile)
    a = systemfunc("rm "+newfile)
  end if
  print("Creating file: "+newfile)
  saveFile = addfile(newfile, "c")
  print("saving date...")
  saveFile->date=date
  print("saving ozone...")
  saveFile->O3_col=O3_col
  print("saving Z0...")
  saveFile->Z0=Z0
  print("Done!")

end


;historyPath="/Volumes/Data/Models/waccm/run/carma-soot/"
;historyRun="Soot5Tg"
;subDir="h0"

;saveMeanFiles(historyRun,historyPath,subDir)

;historyRun="PassiveSoot"

;historyPath="$MODEL_PATH/ccsm/run/"
historyPath="/data5/mmills/"
;historyRun="b40.rcp4_5.2deg.wcm.carma.bc5tgIndPak.001"
;historyRun="b40.rcp4_5.2deg.wcm.001"
;historyRun="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.003"
historyRun="b40.rcp4_5.2deg.wcm.carma.bc0tg.002"
subDir="h0zm"

saveMeanFiles(historyRun,historyPath,subDir)

