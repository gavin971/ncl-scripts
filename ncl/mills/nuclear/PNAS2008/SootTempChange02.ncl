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
  delete(lat)

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
     
  print("reading Soot...")
  Soot=calcGlobMean(history, historyNames, "SOOT01", press, wgt)

  print("saving soot...")
  saveFile->SOOT01=Soot
  
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

undef("saveMeanLOxVars")
procedure saveMeanLOxVars(historyRun:string, historyPattern:string, saveFile:file)

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
 
  ; Loss rates in molec/cm3/s:
  POx_JO2=addfiles_GetVar(history, historyNames, "POx_JO2")
  POx_JO2_av=dim_avg_Wrap(POx_JO2(time|:, lev|:, lat|:, lon|:))
  delete(POx_JO2)
  saveFile->POx_JO2_av=POx_JO2_av
  delete(POx_JO2_av)
  
  LOx_CLO_O=addfiles_GetVar(history, historyNames, "LOx_CLO_O")
  LOx_CLO_O_av=dim_avg_Wrap(LOx_CLO_O(time|:, lev|:, lat|:, lon|:))
  delete(LOx_CLO_O)
  saveFile->LOx_CLO_O_av=LOx_CLO_O_av
  delete(LOx_CLO_O_av)
  
  LOx_HO2_O=addfiles_GetVar(history, historyNames, "LOx_HO2_O")
  LOx_HO2_O_av=dim_avg_Wrap(LOx_HO2_O(time|:, lev|:, lat|:, lon|:))
  delete(LOx_HO2_O)
  saveFile->LOx_HO2_O_av=LOx_HO2_O_av
  delete(LOx_HO2_O_av)

  LOx_HO2_O3=addfiles_GetVar(history, historyNames, "LOx_HO2_O3")
  LOx_HO2_O3_av=dim_avg_Wrap(LOx_HO2_O3(time|:, lev|:, lat|:, lon|:))
  delete(LOx_HO2_O3)
  saveFile->LOx_HO2_O3_av=LOx_HO2_O3_av
  delete(LOx_HO2_O3_av)
  
  LOx_H_O3=addfiles_GetVar(history, historyNames, "LOx_H_O3")
  LOx_H_O3_av=dim_avg_Wrap(LOx_H_O3(time|:, lev|:, lat|:, lon|:))
  delete(LOx_H_O3)
  saveFile->LOx_H_O3_av=LOx_H_O3_av
  delete(LOx_H_O3_av)
  
  LOx_NO2_O=addfiles_GetVar(history, historyNames, "LOx_NO2_O")
  LOx_NO2_O_av=dim_avg_Wrap(LOx_NO2_O(time|:, lev|:, lat|:, lon|:))
  delete(LOx_NO2_O)
  saveFile->LOx_NO2_O_av=LOx_NO2_O_av
  delete(LOx_NO2_O_av)
  
  LOx_O3_O=addfiles_GetVar(history, historyNames, "LOx_O3_O")
  LOx_O3_O_av=dim_avg_Wrap(LOx_O3_O(time|:, lev|:, lat|:, lon|:))
  delete(LOx_O3_O)
  saveFile->LOx_O3_O_av=LOx_O3_O_av
  delete(LOx_O3_O_av)
  
  LOx_OH_O=addfiles_GetVar(history, historyNames, "LOx_OH_O")
  LOx_OH_O_av=dim_avg_Wrap(LOx_OH_O(time|:, lev|:, lat|:, lon|:))
  delete(LOx_OH_O)
  saveFile->LOx_OH_O_av=LOx_OH_O_av
  delete(LOx_OH_O_av)
  
  LOx_O_O=addfiles_GetVar(history, historyNames, "LOx_O_O")
  LOx_O_O_av=dim_avg_Wrap(LOx_O_O(time|:, lev|:, lat|:, lon|:))
  delete(LOx_O_O)
  saveFile->LOx_O_O_av=LOx_O_O_av
  delete(LOx_O_O_av)
  
  O3=addfiles_GetVar(history, historyNames, "O3")
  M=addfiles_GetVar(history, historyNames, "M")
  O3_conc=O3*M
  delete(O3)
  delete(M)
  O3_conc_av=dim_avg_Wrap(O3_conc(time|:, lev|:, lat|:, lon|:))
  delete(O3_conc)
  saveFile->O3_conc_av=O3_conc_av
  delete(O3_conc_av)
  
end

undef("saveMeanFiles")
procedure saveMeanFiles(historyRun:string)

local historyPattern1,historyNames1,history1,historyPattern2,historyNames2, history2,T,Tmean,pi,pm,press,temp,tmean,saveFile1, saveFile2

begin
;  saveFile1 = addfile(historyRun+"/"+historyRun+"_mean1.nc", "c")
;  saveFile2 = addfile(historyRun+"/"+historyRun+"_mean2.nc", "c")
  saveFile3 = addfile(historyRun+"/"+historyRun+"_mean3.nc", "c")
  
;  historyPattern1=historyRun+"/"+historyRun+".cam2.h3.202[56]*.nc"
;  historyPattern2=historyRun+"/"+historyRun+".cam2.h3.2027*.nc"
  historyPattern3=historyRun+"/"+historyRun+".cam2.h3.2028*.nc"
  
;  saveMeanVars(historyRun, historyPattern1, saveFile1)
;  saveMeanVars(historyRun, historyPattern2, saveFile2)  
  saveMeanVars(historyRun, historyPattern3, saveFile3)  
  
;  saveMeanLOxCols(historyRun, historyPattern1, saveFile1)
;  saveMeanLOxCols(historyRun, historyPattern2, saveFile2)
  saveMeanLOxCols(historyRun, historyPattern3, saveFile3)

;  saveMeanLOxVars(historyRun, historyPattern1, saveFile1)
;  saveMeanLOxVars(historyRun, historyPattern2, saveFile2)
end

undef("plotTDiff")
procedure plotTDiff(T_x[*][*]:float, T_c[*][*]:float, Soot[*][*]:float)
local T_diff,T_x,T_c,wks,res,lev
begin
  press=T_x&press
  print("pressure range:"+press(0)+" to "+press(dimsizes(press)-1))
  T_diff=T_x-T_c
  T_diff!0 = "time"
  T_diff!1 = "press"
  T_diff&press = press

  wks = gsn_open_wks("pdf", "TempChange")
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
;  res@lbLabelBarOn        = False         ; turn off color labelbar
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 0.8           ; change aspect ratio of plot
  res@vpHeightF           = 0.55
  res@cnLevelSpacingF     = 10.
  res@gsnRightString      = "(K)"
  res@gsnLeftStringFontHeightF = 24
  res@gsnRightStringFontHeightF = 24
  res@tiYAxisString       = "Pressure (hPa)"
;  res@tiXAxisString       = "Days after May 15"
  res@cnLineLabelPlacementMode= "Computed"
  res@cnLineLabelDensityF = 5 
  
;===================================================================
; Do special labels on Y axis. There will be 6 tickmarks per each
; month, with a month label appearing in the middle.  Here, we are
; treating minor tickmarks like major tickmarks, and vice versa, since
; only the major tickmarks can actually be labeled.
;
  months  = new(8+12+12,integer)
  months(0:7) = (/ 5, 6, 7, 8, 9,10,11,12/)
  months(8:19) =(/ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12/) 
  months(20:31) =(/ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12/) 
  years  = new(8+12+12,integer)
  years(0:7) = 2025
  years(8:19) =2026 
  years(20:31) =2027
  tk_per_mon  =6                                      ; Tickmarks per month.
  nmonths     = dimsizes(months)
  ndays       = nmonths * tk_per_mon

  day_in_mon  = new(nmonths,float)           ; Number of days in each month
  do i=0,nmonths-1
    day_in_mon(i) = days_in_month(years(i),months(i))*1.
  end do

  major_ticks = new(ndays,float)    ; Array to hold major tickmark values.
  imonth      = 0
  days_accum  = 0.
  ii5         = ispan(1,5,1)*1.
  do i=0,ndays-1,6
    major_ticks(i) = days_accum
    major_ticks(i+1:i+5) = major_ticks(i) + ii5*(day_in_mon(imonth)/tk_per_mon)
    days_accum     = days_accum + day_in_mon(imonth)
    imonth = imonth + 1
  end do

  res@tmXBMode           = "Explicit"         ; Define own tick mark labels.
  res@tmXBMinorValues    = major_ticks(3:ndays-1:6)
  res@tmXBValues         = major_ticks
  res@tmXBLabels         = new(ndays,string)
  res@tmXBLabels(:)     = ""

 ; Label mid tickmark with the month label.
  res@tmXBLabels(::6)  = (/"M","J","J","A","S","O","N","D", "J","F","M","A","M","J","J","A","S","O","N","D", "J","F","M","A","M","J","J","A","S","O","N","D"/)

  res@tmXBMinorOutwardLengthF = 0.02    ; Increase length and size of
  res@tmXBMinorLengthF        = 0.02    ; minor tickmarks
  res@tmXBMinorThicknessF     = 2.0
  res@tmXBMajorOutwardLengthF = 0.01
  res@tmXBMajorLengthF        = 0.01
;===================================================================
  
  res@tmYLMode   = "Explicit" 
  res@tmYLValues  = (/    1E-5 ,    3E-5   ,   1E-4  ,    3E-4   , 0.001 , 0.003 , 0.01 ,  0.03 , 0.1 , 0.3 , 1., 3., 10., 30., 100., 300./)
  res@tmYLLabels  = (/"10~S~-5","3x10~S~-5","10~S~-4","3x10~S~-4","0.001","0.003","0.01", "0.03","0.1","0.3","1","3","10","30","100","300"/)
;  res@tmLabelAutoStride = True  
 
;  res@gsnLeftString       = "Soot run, Temperature"
;  plot = gsn_csm_pres_hgt (wks,T_x(press|:, time|:),res)   ; contour the variable
;  res@gsnLeftString       = "Control run, Temperature"
;  plot = gsn_csm_pres_hgt(wks,T_c(press|:, time|:),res)      ; contour the variable
  
  res@gsnLeftString       = "Temperature, Soot - Control, global average"
  res@cnLineLabelInterval = 1 ; label every contour
;  res@cnMonoLevelFlag = True
;  res@cnLevelFlag = "LineAndLabel"
  res@cnLineLabelPlacementMode = "Computed"
;  res@cnLineDashSegLenF = 0.05
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnMinLevelValF=-20
  res@cnMaxLevelValF=100
  res@cnLevelSpacingF     = 10.
  plot = gsn_csm_pres_hgt(wks,T_diff(press|:, time|:),res)      ; contour the variable
  
; Define our own contour levels.  
  res@cnLevelSelectionMode = "ExplicitLevels" 
  res@cnLevels  = (/1,2,3,5,7,10,20,30,50,70,100/) 
;  res@cnLineLabelStrings = (/"1","2","3","5","7","10","20","30","50","70","100"/) 
  Soot=Soot*1e9
  res@gsnLeftString       = "Soot mmr, global average"
  res@gsnRightString      = "(10~S~-9~N~ kg/kg air)"
  plot = gsn_csm_pres_hgt(wks,Soot(press|:, time|:),res)      ; contour the variable
end


undef("plotO3Diff")
procedure plotO3Diff(O3_col_x[*][*]:float, O3_col_c[*][*]:float) 
local O3_col_diff,wks,res,lev
begin 

  lat=O3_col_x&lat
  
  O3_col_diff=O3_col_x-O3_col_c
  O3_col_diff!0 = "time"
  O3_col_diff!1 = "lat"
  O3_col_diff&lat = lat

  O3_col_dpct=100.*O3_col_diff/O3_col_c
  O3_col_dpct!0 = "time"
  O3_col_dpct!1 = "lat"
  O3_col_dpct&lat = lat
  
  wks = gsn_open_wks("pdf", "O3ColumnChange")
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
;  res@lbLabelBarOn        = False         ; turn off color labelbar
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 0.8           ; change aspect ratio of plot
  res@vpHeightF           = 0.55
  res@gsnLeftStringFontHeightF = 24
  res@gsnRightStringFontHeightF = 24
  res@tiYAxisString       = "Latitude"
  res@cnLineLabelPlacementMode= "Computed"
  res@cnLineLabelDensityF = 5 
  
;===================================================================
; Do special labels on Y axis. There will be 6 tickmarks per each
; month, with a month label appearing in the middle.  Here, we are
; treating minor tickmarks like major tickmarks, and vice versa, since
; only the major tickmarks can actually be labeled.
;
  months  = new(8+12+12,integer)
  months(0:7) = (/ 5, 6, 7, 8, 9,10,11,12/)
  months(8:19) =(/ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12/) 
  months(20:31) =(/ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12/) 
  years  = new(8+12+12,integer)
  years(0:7) = 2025
  years(8:19) =2026 
  years(20:31) =2027
  tk_per_mon  =6                                      ; Tickmarks per month.
  nmonths     = dimsizes(months)
  ndays       = nmonths * tk_per_mon

  day_in_mon  = new(nmonths,float)           ; Number of days in each month
  do i=0,nmonths-1
    day_in_mon(i) = days_in_month(years(i),months(i))*1.
  end do

  major_ticks = new(ndays,float)    ; Array to hold major tickmark values.
  imonth      = 0
  days_accum  = 0.
  ii5         = ispan(1,5,1)*1.
  do i=0,ndays-1,6
    major_ticks(i) = days_accum
    major_ticks(i+1:i+5) = major_ticks(i) + ii5*(day_in_mon(imonth)/tk_per_mon)
    days_accum     = days_accum + day_in_mon(imonth)
    imonth = imonth + 1
  end do

  res@tmXBMode           = "Explicit"         ; Define own tick mark labels.
  res@tmXBMinorValues    = major_ticks(3:ndays-1:6)
  res@tmXBValues         = major_ticks
  res@tmXBLabels         = new(ndays,string)
  res@tmXBLabels(:)     = ""

 ; Label mid tickmark with the month label.
  res@tmXBLabels(::6)  = (/"M","J","J","A","S","O","N","D", "J","F","M","A","M","J","J","A","S","O","N","D", "J","F","M","A","M","J","J","A","S","O","N","D"/)

  res@tmXBMinorOutwardLengthF = 0.02    ; Increase length and size of
  res@tmXBMinorLengthF        = 0.02    ; minor tickmarks
  res@tmXBMinorThicknessF     = 2.0
  res@tmXBMajorOutwardLengthF = 0.01
  res@tmXBMajorLengthF        = 0.01
;===================================================================
    
; Define our own contour levels.  
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnMinLevelValF=120
  res@cnMaxLevelValF=510
  res@cnLevelSpacingF     = 30.
  res@gsnRightString      = "O3 Column(DU)"
  res@lbLabelAngleF       = 60
  res@gsnLeftString       = "Control run"
  plot = gsn_csm_lat_time(wks,O3_col_c(lat|:, time|:),res)      ; contour the variable
  res@gsnLeftString       = "Soot run"
  plot = gsn_csm_lat_time(wks,O3_col_x(lat|:, time|:),res)      ; contour the variable
  
;  res@cnLevelSelectionMode = "ManualLevels" 
;  res@cnMinLevelValF=-330
;  res@cnMaxLevelValF=180
  res@cnLevelSelectionMode = "ExplicitLevels" 
  res@cnLevels          = (/-350,-300,-250,-200,-150,-100,-75, -50,-40,-30,-20,-10,-5,-2,0,10,20,40,80,120,160/) 
  res@cnLineThicknesses = (/  1 ,  1 ,  1 ,  1 ,  1 ,  1 , 1 ,  1 , 1 , 1 , 1 , 1 , 1, 1,2, 1, 1, 1, 1, 1 , 1/)
  res@cnMonoLineThickness = False
  res@cnLevelSpacingF     = 30.
  res@gsnLeftString       = "Soot run - Control run"
  plot = gsn_csm_lat_time(wks,O3_col_diff(lat|:, time|:),res)      ; contour the variable

  res@cnLevelSelectionMode = "ExplicitLevels" 
  delete(res@cnLevels)
  res@cnLevels          = (/-60,-50,-40,-30,-20,-10,-5,-2,0,10,20,40,80,120/)
  delete(res@cnLineThicknesses)
  res@cnLineThicknesses = (/ 1 , 1 , 1 , 1 , 1 , 1 , 1, 1,2, 1, 1, 1, 1, 1 /)
  res@cnMonoLineThickness = False
  res@gsnLeftString       = "(Soot - Control)/Control"
  res@gsnRightString      = "O3 Column Change(%)"
  plot = gsn_csm_lat_time(wks,O3_col_dpct(lat|:, time|:),res)      ; contour the variable
end

undef("saveMeans")
procedure saveMeans(historyRun:string, controlRun:string)

begin
  
  saveMeanFiles(historyRun)
  
  saveMeanFiles(controlRun)

end

undef("plotTempO3Change")
procedure plotTempO3Change(historyRun:string, controlRun:string)
local T_x,T_c

begin
  print((/historyRun/))
;  saveMeanVars(historyRun,True)
;  saveMeanLOxFiles(historyRun)
  historyPattern = historyRun+"/"+historyRun+"_mean*.nc"
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  T_x=addfiles_GetVar(history, historyNames, "T")
  Soot=addfiles_GetVar(history, historyNames, "SOOT01")
  O3_col_x=addfiles_GetVar(history, historyNames, "O3_col")
  O3_col_x=O3_col_x/2.687e16 ; Convert to Dobson Units
  O3_col_x@units="Dobson Units"
 
;  T_x=File1->T
;  Soot=File1->SOOT01

  print((/controlRun/))
;  saveMeanVars(controlRun,False)
;  saveMeanLOxFiles(controlRun)
  controlPattern = controlRun+"/"+controlRun+"_mean*.nc"
  controlNames = systemfunc("ls " + controlPattern)
  control = addfiles(controlNames, "r")
  ListSetType(control, "cat")
  T_c=addfiles_GetVar(control, controlNames, "T")
  O3_col_c=addfiles_GetVar(control, controlNames, "O3_col")
  O3_col_c=O3_col_c/2.687e16 ; Convert to Dobson Units
  O3_col_c@units="Dobson Units"
  
;  File2 = addfile(controlRun+"/"+controlRun+"Tmean.nc", "r")
;  T_c=File2->T

  print("plotting Temperature...")
  plotTDiff(T_x,T_c,Soot)
  print("plotting ozone column...")
  plotO3Diff(O3_col_x, O3_col_c)
    
end

undef("LOxLatMean")
function LOxLatMean(Rate[*][*]:float,j:integer)

local pi,wgt
begin
  pi=acos(0.0)*2.0
  if (j.eq.-1) then
    Rate!1 = "lat"
    lat=Rate&lat
    wgt=dble2flt(cos(lat*pi/180.))
    LOx=dim_avg_wgt_Wrap(Rate,wgt,0)
    delete(wgt)
  else
    LOx=Rate(:,j)
  end if
  delete(pi)
return(LOx)
end

undef("plotLOx")
procedure plotLOx(historyRun:string, controlRun:string)
local T_x,T_c

begin
  print((/historyRun/))
  historyPattern = historyRun+"/"+historyRun+"_mean*.nc"
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  
  R_BrO_BrO_col_x = addfiles_GetVar(history, historyNames, "R_BrO_BrO_col")
  R_BrO_ClOa_col_x = addfiles_GetVar(history, historyNames, "R_BrO_ClOa_col")
  R_BrO_ClOb_col_x = addfiles_GetVar(history, historyNames, "R_BrO_ClOb_col")
  R_BrO_O_col_x = addfiles_GetVar(history, historyNames, "R_BrO_O_col")
  R_CH3O2_NO_col_x = addfiles_GetVar(history, historyNames, "R_CH3O2_NO_col")
  R_ClO_O_col_x = addfiles_GetVar(history, historyNames, "R_CLO_O_col")
  R_HO2_O_col_x = addfiles_GetVar(history, historyNames, "R_HO2_O_col")
  R_HO2_O3_col_x = addfiles_GetVar(history, historyNames, "R_HO2_O3_col")
  R_H_O3_col_x = addfiles_GetVar(history, historyNames, "R_H_O3_col")
  R_NO2_O_col_x = addfiles_GetVar(history, historyNames, "R_NO2_O_col")
  R_NO_HO2_col_x = addfiles_GetVar(history, historyNames, "R_NO_HO2_col")
  R_NO_N_col_x = addfiles_GetVar(history, historyNames, "R_NO_N_col")
  R_N_O2_col_x = addfiles_GetVar(history, historyNames, "R_N_O2_col")
  R_O1D_H2O_col_x = addfiles_GetVar(history, historyNames, "R_O1D_H2O_col")
  R_O2_O_col_x = addfiles_GetVar(history, historyNames, "R_O2_O_col")
  R_O3_O_col_x = addfiles_GetVar(history, historyNames, "R_O3_O_col")
  R_OH_O_col_x = addfiles_GetVar(history, historyNames, "R_OH_O_col")
  R_O_O_col_x = addfiles_GetVar(history, historyNames, "R_O_O_col")

  J_Cl2O2_col_x = addfiles_GetVar(history, historyNames, "J_Cl2O2_col")
  J_HOBr_col_x = addfiles_GetVar(history, historyNames, "J_HOBr_col")
  J_HOCl_col_x = addfiles_GetVar(history, historyNames, "J_HOCl_col")
  J_NO_col_x = addfiles_GetVar(history, historyNames, "J_NO_col")
  J_NO3_b_col_x = addfiles_GetVar(history, historyNames, "J_NO3_b_col")
  J_O2_a_col_x = addfiles_GetVar(history, historyNames, "J_O2_a_col")
  J_O2_b_col_x = addfiles_GetVar(history, historyNames, "J_O2_b_col")

  O3_col_x=addfiles_GetVar(history, historyNames, "O3_col")
  O3_col_x=O3_col_x/2.687e16 ; Convert to Dobson Units
  O3_col_x@units="Dobson Units"

  print((/controlRun/))
;  saveMeanVars(controlRun,False)
;  saveMeanLOxFiles(controlRun)
  controlPattern = controlRun+"/"+controlRun+"_mean*.nc"
  controlNames = systemfunc("ls " + controlPattern)
  control = addfiles(controlNames, "r")
  ListSetType(control, "cat")
  R_BrO_BrO_col_c = addfiles_GetVar(control, controlNames, "R_BrO_BrO_col")
  R_BrO_ClOa_col_c = addfiles_GetVar(control, controlNames, "R_BrO_ClOa_col")
  R_BrO_ClOb_col_c = addfiles_GetVar(control, controlNames, "R_BrO_ClOb_col")
  R_BrO_O_col_c = addfiles_GetVar(control, controlNames, "R_BrO_O_col")
  R_CH3O2_NO_col_c = addfiles_GetVar(control, controlNames, "R_CH3O2_NO_col")
  R_ClO_O_col_c = addfiles_GetVar(control, controlNames, "R_CLO_O_col")
  R_HO2_O_col_c = addfiles_GetVar(control, controlNames, "R_HO2_O_col")
  R_HO2_O3_col_c = addfiles_GetVar(control, controlNames, "R_HO2_O3_col")
  R_H_O3_col_c = addfiles_GetVar(control, controlNames, "R_H_O3_col")
  R_NO2_O_col_c = addfiles_GetVar(control, controlNames, "R_NO2_O_col")
  R_NO_HO2_col_c = addfiles_GetVar(control, controlNames, "R_NO_HO2_col")
  R_NO_N_col_c = addfiles_GetVar(control, controlNames, "R_NO_N_col")
  R_N_O2_col_c = addfiles_GetVar(control, controlNames, "R_N_O2_col")
  R_O1D_H2O_col_c = addfiles_GetVar(control, controlNames, "R_O1D_H2O_col")
  R_O2_O_col_c = addfiles_GetVar(control, controlNames, "R_O2_O_col")
  R_O3_O_col_c = addfiles_GetVar(control, controlNames, "R_O3_O_col")
  R_OH_O_col_c = addfiles_GetVar(control, controlNames, "R_OH_O_col")
  R_O_O_col_c = addfiles_GetVar(control, controlNames, "R_O_O_col")

  J_Cl2O2_col_c = addfiles_GetVar(control, controlNames, "J_Cl2O2_col")
  J_HOBr_col_c = addfiles_GetVar(control, controlNames, "J_HOBr_col")
  J_HOCl_col_c = addfiles_GetVar(control, controlNames, "J_HOCl_col")
  J_NO_col_c = addfiles_GetVar(control, controlNames, "J_NO_col")
  J_NO3_b_col_c = addfiles_GetVar(control, controlNames, "J_NO3_b_col")
  J_O2_a_col_c = addfiles_GetVar(control, controlNames, "J_O2_a_col")
  J_O2_b_col_c = addfiles_GetVar(control, controlNames, "J_O2_b_col")

  O3_col_c=addfiles_GetVar(control, controlNames, "O3_col")
  O3_col_c=O3_col_c/2.687e16 ; Convert to Dobson Units
  O3_col_c@units="Dobson Units"
    
  time=R_O3_O_col_x&time
  lat=R_O3_O_col_x&lat
 
  ; Create a panel plot of the loss rates.
  print("")
  j=-1
  res = True
  resP = True
  if (j.eq.-1) then
    print(" Plotting global mean...")
    wks = gsn_open_wks("pdf", "LOx_GlobAvg")
    resP@txString = "Global Mean"
  else
    print(" Plotting at latitude "+lat(j)+"...")
    wks = gsn_open_wks("pdf", "LOx_"+lat(j))
    resP@txString = "Latitude: "+lat(j)
  end if
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
    
;  res@gsnMaximize           = True
  res@vpWidthF              = 1.0
  res@gsnDraw = False
  res@gsnFrame = False
  
  res@trYMinF = 0
;  res@trYMaxF = 5.0
  
;  res@tmYLFormat = "0*+^e"
;  res@tmYLPrecision = 3
  res@pmLegendDisplayMode = "Always"
  res@pmLegendWidthF = 0.4
;------Page 1--------  

  nPanel=3
  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)
  
  data= new((/2,dimsizes(time)/),float)
  
  res@gsnRightString = "DU"
  res@gsnLeftString = "O3 column"
  data= new((/2,dimsizes(time)/),float)
  data(0,:)=LOxLatMean(O3_col_x,j)
  data(1,:)=LOxLatMean(O3_col_c,j)
  res@xyExplicitLegendLabels=(/"Soot","Control"/)
  plot(0) = gsn_csm_xy(wks, time, data, res)  
  
  delete(data)
  res@gsnLeftString = "Control: Ox Column Production & Loss"
  res@gsnRightString = "molec/cm~S~2~N~/s"
  data= new((/5,dimsizes(time)/),float)
  data(0,:)=2*LOxLatMean(R_O3_O_col_c,j)
  data(1,:)=2*LOxLatMean(R_ClO_O_col_c,j)
  data(2,:)=2*LOxLatMean(R_O_O_col_c,j)
  data(3,:)=2*LOxLatMean(R_NO2_O_col_c,j)
  Rate=J_O2_a_col_c+J_O2_b_col_c
  Rate!1="lat"
  Rate&lat=lat
  data(4,:)=LOxLatMean(Rate,j)
  delete(Rate)
  res@xyLineThicknessF = 4.0
  res@xyLineColors      = (/"blue","red","green","black","orange"/); change line color
  delete(res@xyExplicitLegendLabels)
  res@xyExplicitLegendLabels=(/"2 x O~B~3~N~+O","2 x ClO+O","2 x O+O","2 x NO~B~2~N~+O","O~B~2~N~+hv"/)
  plot(1) = gsn_csm_xy(wks, time, data, res)  
  
  delete(data)
  res@pmLegendDisplayMode = "NoCreate"
  res@gsnLeftString = "Soot: Ox Column Production & Loss"
  res@gsnRightString = "molec/cm~S~2~N~/s"
  data= new((/5,dimsizes(time)/),float)
  data(0,:)=2*LOxLatMean(R_O3_O_col_x,j)
  data(1,:)=2*LOxLatMean(R_ClO_O_col_x,j)
  data(2,:)=2*LOxLatMean(R_O_O_col_x,j)
  data(3,:)=2*LOxLatMean(R_NO2_O_col_x,j)
  Rate=J_O2_a_col_x+J_O2_b_col_x
  Rate!1="lat"
  Rate&lat=lat
  data(4,:)=LOxLatMean(Rate,j)
  delete(Rate)
  res@xyLineThicknessF = 4.0
  res@xyLineColors      = (/"blue","red","green","black","orange"/); change line color
  plot(2) = gsn_csm_xy(wks, time, data, res)  

  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(plot)

;------Page 2--------  
  res@pmLegendDisplayMode = "Always"
  
  nPanel=6
  res@vpHeightF             = 1.0 / nPanel              ; change aspect ratio
  plot = new(nPanel, graphic)
  
  res@xyLineThicknessF = 4.0
  
  delete(res@xyLineColors)
  res@xyLineColors      = (/"blue","red"/)          ; change line color
  
  delete(data)
  data= new((/2,dimsizes(time)/),float)
    
  res@gsnRightString = "molec/cm~S~2~N~/s"
  delete(res@xyExplicitLegendLabels)
  res@xyExplicitLegendLabels=(/"Soot","Control"/)
    
  res@gsnLeftString = "O2+hv rate"
  Rate=J_O2_a_col_x+J_O2_b_col_x
  Rate!1="lat"
  Rate&lat=lat
  data(0,:)=LOxLatMean(Rate,j)
  Rate=J_O2_a_col_c+J_O2_b_col_c
  data(1,:)=LOxLatMean(Rate,j)
  delete(Rate)
  plot(0) = gsn_csm_xy(wks, time, data, res)
  
  res@pmLegendDisplayMode = "NoCreate"
  res@gsnLeftString = "2 x O3+O rate"
  data(0,:)=2*LOxLatMean(R_O3_O_col_x,j)
  data(1,:)=2*LOxLatMean(R_O3_O_col_c,j)
  plot(1) = gsn_csm_xy(wks, time, data, res)

  res@gsnLeftString = "2 x ClO+O rate"
  data(0,:)=2*LOxLatMean(R_ClO_O_col_x,j)
  data(1,:)=2*LOxLatMean(R_ClO_O_col_c,j)
  plot(2) = gsn_csm_xy(wks, time, data, res)
  
  res@gsnLeftString = "2 x O+O rate"
  data(0,:)=2*LOxLatMean(R_O_O_col_x,j)
  data(1,:)=2*LOxLatMean(R_O_O_col_c,j)
  plot(3) = gsn_csm_xy(wks, time, data, res)
  
  res@gsnLeftString = "2 x NO2+O rate"
  data(0,:)=2*LOxLatMean(R_NO2_O_col_x,j)
  data(1,:)=2*LOxLatMean(R_NO2_O_col_c,j)
  plot(4) = gsn_csm_xy(wks, time, data, res)
  
  res@gsnLeftString = "O2+O rate"
  data(0,:)=LOxLatMean(R_O2_O_col_x,j)
  data(1,:)=LOxLatMean(R_O2_O_col_c,j)
  plot(5) = gsn_csm_xy(wks, time, data, res)
  
  gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
  delete(wks)
 
    
end
