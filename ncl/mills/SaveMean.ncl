load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "tests/getHistoryPath.ncl"

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

undef("calcZdel")
function calcZdel(history:list, historyNames[*]:string)

local z3,zdel
begin
 
  print("Reading geopotential heights...")
  z3 = addfiles_GetVar(history, historyNames, "Z3")
  nlev=dimsizes(z3&lev)
  
  print("Calculating zdel...")
  zdel = z3
  zdel(:,0:nlev-2,:,:)=z3(:,0:nlev-2,:,:)-z3(:,1:nlev-1,:,:)
    
  delete(z3)
  delete(nlev)
  
  return(zdel)
    
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
  columnAv=dim_avg_Wrap(column(time|:, lat|:, lon|:))
  delete(column)
  columnAv!0 = "time"
  columnAv!1 = "lat"
  columnAv&lat = lat
  delete(lat)

  return(columnAv)
end

undef("calcCol_conc")
function calcCol_conc(history:list, historyNames[*]:string, var:string, zdel:float)

local conc,levcol,column,columnAv
begin
  ; Some constants we will need:
  g       = 980.6               ; cm s^-2 
  Mair    = 28.964/6.022e23     ; g/molec
  pi=acos(0.0)*2.0
  print("calculating column")

  ; volume mixing ratios:
  print("Reading vmr...")
  conc=addfiles_GetVar(history, historyNames, var)
  levcol=conc*zdel
  lat=conc&lat
  delete(conc)
  levcol!0 = "time"
  levcol!1 = "lev"
  levcol!2 = "lat"
  levcol!3 = "lon"

  print("Calculating column ...")
  column = dim_sum_Wrap(levcol(time|:, lat|:, lon|:, lev|:))
  delete(levcol)
  columnAv=dim_avg_Wrap(column(time|:, lat|:, lon|:))
  delete(column)
  columnAv!0 = "time"
  columnAv!1 = "lat"
  columnAv&lat = lat
;  delete(lat)

  return(columnAv)
end

undef("calcGlobMean")
function calcGlobMean(history:list, historyNames[*]:string, var:string, press[*]:float, wgt[*]:float)

local temp,tmean
begin

  print("reading "+var)
  temp=addfiles_GetVar(history, historyNames, var)
  print(dimsizes(temp))
  
  print("averaging zonally")
  tmean=dim_avg_Wrap(temp(time|:, lev|:, lat|:, lon|:))
  
  delete(temp)
  
  print("averaging meridionally")
  temp=dim_avg_wgt_Wrap(tmean(time|:, lev|:, lat|:),wgt,0)
  temp!0 = "time"
  temp!1 = "press"
  temp&press = press
  
  return(temp)

end

undef("saveMeanVars")
procedure saveMeanVars(historyRun:string, historyPattern:string, saveFile:file)

local historyPattern,historyNames,history,T,Tmean,pi,pm,press,temp,tmean

begin

  pi=acos(0.0)*2.0
   
  ; Open the set of history files.
  ;
  ; NOTE: We want to work with the entire set of history files, 
  ; so we can do a time evolution.
  print("  Reading data ...")
  
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  print("historyPattern: "+historyPattern)
  print("historyNames:"+dimsizes(historyNames))
  print(historyNames)
    
  
  ; Read in the fields we will need.
 ; time = addfiles_GetVar(history, historyNames, "time")
  lev = history[0]->lev
  lat = history[0]->lat
  lon = history[0]->lon
  wgt=dble2flt(cos(lat*pi/180.))
  
  ; Ozone:
  print("reading ozone...")
  press=calcPress(history, historyNames)
  pdel=calcPdel(history, historyNames)
  O3_col=calcCol_vmr(history, historyNames, "O3", pdel)
  O3_col@units="molec/cm2"
  delete(pdel)
  
  print("saving ozone...")
  saveFile->O3_col=O3_col
  
  ; Temperature:
  print("reading temperature...")
  temp=calcGlobMean(history, historyNames, "T", press, wgt)
     
  print("saving temperature...")
    
  saveFile->T=temp
     
  print("reading SAD_SULF...")
  SAD_SULF=calcGlobMean(history, historyNames, "SAD_SULF", press, wgt)

  print("saving SAD_SULF...")
  saveFile->SAD_SULF=SAD_SULF
  
  print("done saving SAD_SULF...")
end

undef("saveMeanLOxCols")
procedure saveMeanLOxCols(historyRun:string, historyPattern:string, saveFile:file)

local historyNames,history,saveFile

begin

  pi=acos(0.0)*2.0
   
  ; Open the set of history files.
  ;
  ; NOTE: We want to work with the entire set of history files, 
  ; so we can do a time evolution.
  print("  Reading data ...")
  
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  print("historyNames:"+dimsizes(historyNames))
 ; print(historyNames1)
      
 ; Read in the fields we will need.
  lev = history[0]->lev
  lat = history[0]->lat
  lon = history[0]->lon
 
  print("calculating zdel")
  zdel=calcZdel(history, historyNames)

  ; Loss rates in molec/cm3/s:
  print("J_Cl2O2")
  J_Cl2O2_col   =calcCol_conc(history, historyNames, "J_Cl2O2", zdel)
  J_Cl2O2_col@units="molec/cm2/s"
  saveFile->J_Cl2O2_col=J_Cl2O2_col
  delete(J_Cl2O2_col)
  
  print("J_HOBr")
  J_HOBr_col =calcCol_conc(history, historyNames, "J_HOBr", zdel)
  saveFile->J_HOBr_col=J_HOBr_col
  delete(J_HOBr_col)
  
  print("J_HOCl")
  J_HOCl_col =calcCol_conc(history, historyNames, "J_HOCl", zdel)
  saveFile->J_HOCl_col=J_HOCl_col
  delete(J_HOCl_col)
    
  print("J_NO")
  J_NO_col=calcCol_conc(history, historyNames, "J_NO", zdel)
  saveFile->J_NO_col=J_NO_col
  delete(J_NO_col)
  
  print("J_NO3_b")
  J_NO3_b_col  =calcCol_conc(history, historyNames, "J_NO3_b", zdel)
  saveFile->J_NO3_b_col=J_NO3_b_col
  delete(J_NO3_b_col)
  
  print("J_O2_a")
  J_O2_a_col =calcCol_conc(history, historyNames, "J_O2_a", zdel)
  saveFile->J_O2_a_col=J_O2_a_col
  delete(J_O2_a_col)
  
  print("J_O2_b")
  J_O2_b_col  =calcCol_conc(history, historyNames, "J_O2_b", zdel)
  saveFile->J_O2_b_col=J_O2_b_col
  delete(J_O2_b_col)
  
  print("R_BrO_BrO")
  R_BrO_BrO_col  =calcCol_conc(history, historyNames, "R_BrO_BrO", zdel)
  saveFile->R_BrO_BrO_col=R_BrO_BrO_col
  delete(R_BrO_BrO_col)
  
  print("R_BrO_ClOa")
  R_BrO_ClOa_col   =calcCol_conc(history, historyNames, "R_BrO_ClOa", zdel)
  saveFile->R_BrO_ClOa_col=R_BrO_ClOa_col
  delete(R_BrO_ClOa_col)

  print("R_BrO_ClOb")
  R_BrO_ClOb_col   =calcCol_conc(history, historyNames, "R_BrO_ClOb", zdel)
  saveFile->R_BrO_ClOb_col=R_BrO_ClOb_col
  delete(R_BrO_ClOb_col)

  print("R_BrO_O")
  R_BrO_O_col   =calcCol_conc(history, historyNames, "R_BrO_O", zdel)
  saveFile->R_BrO_O_col=R_BrO_O_col
  delete(R_BrO_O_col)

  print("R_CH3O2_NO")
  R_CH3O2_NO_col   =calcCol_conc(history, historyNames, "R_CH3O2_NO", zdel)
  saveFile->R_CH3O2_NO_col=R_CH3O2_NO_col
  delete(R_CH3O2_NO_col)

  print("R_CLO_O")
  R_CLO_O_col   =calcCol_conc(history, historyNames, "R_CLO_O", zdel)
  saveFile->R_CLO_O_col=R_CLO_O_col
  delete(R_CLO_O_col)

  print("R_HO2_O")
  R_HO2_O_col   =calcCol_conc(history, historyNames, "R_HO2_O", zdel)
  saveFile->R_HO2_O_col=R_HO2_O_col
  delete(R_HO2_O_col)

  print("R_HO2_O3")
  R_HO2_O3_col   =calcCol_conc(history, historyNames, "R_HO2_O3", zdel)
  saveFile->R_HO2_O3_col=R_HO2_O3_col
  delete(R_HO2_O3_col)

  print("R_H_O3")
  R_H_O3_col   =calcCol_conc(history, historyNames, "R_H_O3", zdel)
  saveFile->R_H_O3_col=R_H_O3_col
  delete(R_H_O3_col)

  print("R_NO2_O")
  R_NO2_O_col   =calcCol_conc(history, historyNames, "R_NO2_O", zdel)
  saveFile->R_NO2_O_col=R_NO2_O_col
  delete(R_NO2_O_col)

  print("R_NO_HO2")
  R_NO_HO2_col   =calcCol_conc(history, historyNames, "R_NO_HO2", zdel)
  saveFile->R_NO_HO2_col=R_NO_HO2_col
  delete(R_NO_HO2_col)

  print("R_NO_ON")  ; note typo in model output
  R_NO_N_col   =calcCol_conc(history, historyNames, "R_NO_ON", zdel)
  saveFile->R_NO_N_col=R_NO_N_col
  delete(R_NO_N_col)

  print("R_N_O2")
  R_N_O2_col   =calcCol_conc(history, historyNames, "R_N_O2", zdel)
  saveFile->R_N_O2_col=R_N_O2_col
  delete(R_N_O2_col)

  print("R_O1D_H2O")
  R_O1D_H2O_col   =calcCol_conc(history, historyNames, "R_O1D_H2O", zdel)
  saveFile->R_O1D_H2O_col=R_O1D_H2O_col
  delete(R_O1D_H2O_col)

  print("R_O2_O")
  R_O2_O_col   =calcCol_conc(history, historyNames, "R_O2_O", zdel)
  saveFile->R_O2_O_col=R_O2_O_col
  delete(R_O2_O_col)

  print("R_O3_O")
  R_O3_O_col   =calcCol_conc(history, historyNames, "R_O3_O", zdel)
  saveFile->R_O3_O_col=R_O3_O_col
  delete(R_O3_O_col)

  print("R_OH_O")
  R_OH_O_col   =calcCol_conc(history, historyNames, "R_OH_O", zdel)
  saveFile->R_OH_O_col=R_OH_O_col
  delete(R_OH_O_col)

  print("R_O_O")
  R_O_O_col   =calcCol_conc(history, historyNames, "R_O_O", zdel)
  saveFile->R_O_O_col=R_O_O_col
  delete(R_O_O_col)

  delete(zdel)
  
end

undef("saveMeanFiles")
procedure saveMeanFiles(historyRun:string, pattern:string, subDir:string)

local historyPattern1,historyNames1,history1,historyPattern2,historyNames2, \
  history2,T,Tmean,pi,pm,press,temp,tmean,saveFile1,saveFile2

begin

  historyPath = getHistoryPath(historyRun)

;  saveFile1 = addfile(historyRun+"/"+historyRun+".mean.h3.1995.nc", "c")
;  saveFile2 = addfile(historyRun+"/"+historyRun+".mean.h0.1996.nc", "c")
;  saveFile3 = addfile(historyRun+"/"+historyRun+".mean.h3.1997.nc", "c")
  saveFile = addfile(historyPath+historyRun+".mean."+pattern+".nc", "c")
  
;  historyPattern1=historyRun+"/"+historyRun+".cam2.h3.1995*.nc"
;  historyPattern2=historyRun+"/"+historyRun+".cam2.h0.1996*.nc"
;  historyPattern3=historyRun+"/"+historyRun+".cam2.h3.1997*.nc"
  historyPattern=historyPath+subDir+"/*.cam2."+pattern+"*.nc"
  
;  saveMeanVars(historyRun, historyPattern1, saveFile1)
;  saveMeanVars(historyRun, historyPattern2, saveFile2)  
;  saveMeanVars(historyRun, historyPattern3, saveFile3)  
  saveMeanVars(historyRun, historyPattern, saveFile)
  
;  saveMeanLOxCols(historyRun, historyPattern1, saveFile1)
;  saveMeanLOxCols(historyRun, historyPattern2, saveFile2)
;  saveMeanLOxCols(historyRun, historyPattern3, saveFile3)
  saveMeanLOxCols(historyRun, historyPattern, saveFile)
  
end

undef("saveMeans")
procedure saveMeans(historyRun:string, controlRun:string, \
                     pattern:string, subDir:string)

begin

  print((/historyRun/))  
  saveMeanFiles(historyRun, pattern, subDir)
  
  print((/controlRun/))  
  saveMeanFiles(controlRun, pattern, subDir)

end

;Main

;saveMeans("geoengn4","carmasadn3")
;saveMeans("geo2","nongeo2")
saveMeanFiles("nogeoPurSpin1","h3.1995","daily")
;saveMeans("geoPurSpin1","nogeoPurSpin1","h3.1995","daily")
