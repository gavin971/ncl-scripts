load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("calcCol")
function calcCol(pdel[*][*][*][*]:float,vmr[*][*][*][*]:float)

local g,Mair,partPress,column,columnAv
begin
  ; Some constants we will need:
  g       = 980.6               ; cm s^-2 
  Mair    = 28.964/6.022e23     ; g/molec

  ; dP=rho*g*dh=[M]*Mair*g*dh
  ; [M]*dh=dP/(Mair*g)
  partPress=pdel*vmr*10./(Mair*g)  ; pdel(Pa)*10(dynes/cm2/Pa)
  partPress!0 = "time"
  partPress!1 = "lev"
  partPress&lev = vmr&lev
  partPress!2 = "lat"
  partPress&lat = vmr&lat
  partPress!3 = "lon"
  partPress&lon = vmr&lon
  
;  print(pdel(0,:,0,0))
;  print(partPress(0,:,0,0))
  print("Calculating column ...")
  column = dim_sum_Wrap(partPress(time|:, lat|:, lon|:, lev|:))
  columnAv=dim_avg_Wrap(column(time|:, lat|:, lon|:))
  
  return(columnAv)
end

undef("calcO3Col")
function calcO3Col(historyRun:string)

local ps,p0,hyai,hybi,pm,pdel,O3_x,O3_x_col,saveFile

begin

  ; Open the set of history files.
  ;
  ; NOTE: We want to work with the entire set of history files, 
  ; so we can do a time evolution.
  print("  Reading data ...")
  
  historyPattern=historyRun+"/"+historyRun+".cam2.h0.20*.nc"

  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
    
  ; Read in the fields we will need.
  ps = addfiles_GetVar(history, historyNames, "PS")
  p0 = history[0]->P0
  hyai = history[0]->hyai
  hybi = history[0]->hybi
  time = addfiles_GetVar(history, historyNames, "time")
  lev = history[0]->lev
  lat = history[0]->lat
  lon = history[0]->lon
  
  pm = pres_hybrid_ccm(ps,p0,hyai,hybi)
  
  delete(hyai)
  delete(hybi)
;  print(pm(0,:,0,0))
  
  pdel=pm(:,1:dimsizes(lev),:,:)-pm(:,0:dimsizes(lev)-1,:,:)
  
  delete(pm)
  
;  print(pdel(0,:,0,0))
 
  ; Air density:
;  M=addfiles_GetVar(history, historyNames, "M")

  ; Loss rates in molec/cm3/s:
;  POx_JO2_x=addfiles_GetVar(history, historyNames, "POx_JO2")
;  LOx_CLO_O_x=addfiles_GetVar(history, historyNames, "LOx_CLO_O")
;  LOx_HO2_O_x=addfiles_GetVar(history, historyNames, "LOx_HO2_O")
;  LOx_HO2_O3_x=addfiles_GetVar(history, historyNames, "LOx_HO2_O3")
;  LOx_H_O3_x=addfiles_GetVar(history, historyNames, "LOx_H_O3")
;  LOx_NO2_O_x=addfiles_GetVar(history, historyNames, "LOx_NO2_O")
;  LOx_O3_O_x=addfiles_GetVar(history, historyNames, "LOx_O3_O")
;  LOx_OH_O_x=addfiles_GetVar(history, historyNames, "LOx_OH_O")
;  LOx_O_O_x=addfiles_GetVar(history, historyNames, "LOx_O_O")
;
;  POx_JO2_c=addfiles_GetVar(control, controlNames, "POx_JO2")
;  LOx_CLO_O_c=addfiles_GetVar(control, controlNames, "LOx_CLO_O")
;  LOx_HO2_O_c=addfiles_GetVar(control, controlNames, "LOx_HO2_O")
;  LOx_HO2_O3_c=addfiles_GetVar(control, controlNames, "LOx_HO2_O3")
;  LOx_H_O3_c=addfiles_GetVar(control, controlNames, "LOx_H_O3")
;  LOx_NO2_O_c=addfiles_GetVar(control, controlNames, "LOx_NO2_O")
;  LOx_O3_O_c=addfiles_GetVar(control, controlNames, "LOx_O3_O")
;  LOx_OH_O_c=addfiles_GetVar(control, controlNames, "LOx_OH_O")
;  LOx_O_O_c=addfiles_GetVar(control, controlNames, "LOx_O_O")

  ; Ozone mixing ratios:
  O3_x=addfiles_GetVar(history, historyNames, "O3")

  print("calculating O3_x column")
  O3_x_col=calcCol(pdel,O3_x)
  print("max (molec/cm2):")
  print(max((/O3_x_col/)))
  print("min (molec/cm2):")
  print(min((/O3_x_col/)))
  O3_x_col=O3_x_col/2.687e16 ; Convert to Dobson Units
  print("max (DU):")
  print(max((/O3_x_col/)))
  print("min (DU):")
  print(min((/O3_x_col/)))
  O3_x_col&units="Dobson Units"
  
  saveFile = addfile(historyRun+"/"+historyRun+"O3Col.nc", "c")
  saveFile->O3_col=O3_x_col
  
  return(O3_x_col)
end

undef("")
procedure plotOxDiff(O3_x_col[*][*]:float,O3_c_col[*][*]:float)
local O3_diff_col,O3_x_col,O3_c_col,wks,res
begin
  O3_diff_col=O3_x_col-O3_c_col
  O3_diff_col!0 = "time"
  O3_diff_col!1 = "lat"
  O3_diff_col&lat = O3_x_col&lat

  wks = gsn_open_wks("pdf", "O3Column")
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
  res@lbLabelBarOn        = False         ; turn off color labelbar
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 0.8           ; change aspect ratio of plot
  res@vpHeightF           = 0.55
  res@cnLevelSpacingF     = 20.
  res@gsnRightString      = "Dobson Units"
  res@tiYAxisString       = "Latitude"
  res@tiXAxisString       = "Days after May 15"
  
  res@gsnLeftString       = "Soot run, Column O3"
  plot = gsn_csm_contour(wks,O3_x_col(lat|:, time|:),res)      ; contour the variable
  res@gsnLeftString       = "Control run, Column O3"
  plot = gsn_csm_contour(wks,O3_c_col(lat|:, time|:),res)      ; contour the variable
  res@gsnLeftString       = "Soot - Control, Column O3"
  plot = gsn_csm_contour(wks,O3_diff_col(lat|:, time|:),res)      ; contour the variable

end

undef("plotSootOxLoss")
procedure plotSootOxLoss(historyRun:string, controlRun:string)
local O3_diff_col,O3_x_col,O3_c_col

begin
  print((/historyRun/))
;  O3_x_col=calcO3Col(historyRun)
  File1 = addfile(historyRun+"/"+historyRun+"O3Col.nc", "r")
  O3_x_col=File1->O3_col

  print((/controlRun/))
;  O3_c_col=calcO3Col(controlRun)
  File2 = addfile(controlRun+"/"+controlRun+"O3Col.nc", "r")
  O3_c_col=File2->O3_col

  print("plotting O3 column...")
  plotOxDiff(O3_x_col,O3_c_col)
  
end
