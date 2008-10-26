load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("histVarExtract")
procedure histVarExtract(historyPattern:string, saveFile:file)

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
      
 ; Read in the fields we will need.
  print("saving P0, lev, ilev, lat, lon, hyai, hybi, hyam, hybm...")
  saveFile->P0 = history[0]->P0
  saveFile->lev = history[0]->lev
  saveFile->ilev = history[0]->ilev
  saveFile->lat = history[0]->lat
  saveFile->lon = history[0]->lon
  saveFile->hyai = history[0]->hyai
  saveFile->hybi = history[0]->hybi
  saveFile->hyam = history[0]->hyam
  saveFile->hybm = history[0]->hybm
 
  print("saving time...")
  time=addfiles_GetVar(history, historyNames, "time")
  saveFile->time=time
  print("saving date...")
  date=addfiles_GetVar(history, historyNames, "date")
  saveFile->date=date
  print("saving NOX...")
  NOX =addfiles_GetVar(history, historyNames, "NOX")
  saveFile->NOX = NOX
  print("saving NOY...")
  NOY =addfiles_GetVar(history, historyNames, "NOY")
  saveFile->NOY = NOY
  print("saving O3...")
  O3  =addfiles_GetVar(history, historyNames, "O3")
  saveFile->O3  = O3
  print("saving PS...")
  PS  =addfiles_GetVar(history, historyNames, "PS")
  saveFile->PS  =PS
  print("saving PSL...")
  PSL =addfiles_GetVar(history, historyNames, "PSL")
  saveFile->PSL =PSL
  print("saving T...")
  T   =addfiles_GetVar(history, historyNames, "T")
  saveFile->T   =T
  print("saving U...")
  U   =addfiles_GetVar(history, historyNames, "U")
  saveFile->U   =U
  print("saving V...")
  V   =addfiles_GetVar(history, historyNames, "V")
  saveFile->V   =V
  print("saving OMEGA...")
  OMEGA   =addfiles_GetVar(history, historyNames, "OMEGA")
  saveFile->OMEGA   =OMEGA
  print("saving Z3...")
  Z3  =addfiles_GetVar(history, historyNames, "Z3")
  saveFile->Z3  =Z3
  print("saving QSUM...")
  QSUM  =addfiles_GetVar(history, historyNames, "QSUM")
  saveFile->QSUM  =QSUM
;  print("saving TS...")
;  TS  =addfiles_GetVar(history, historyNames, "TS")
;  saveFile->TS  =TS

end

;Main

caseid="noaur_eco"
htype="h3"
path="run/"+caseid+"/"

do year=2027,2040
  historyPattern=path+"*cam2."+htype+"."+year+"*.nc"
;  historyPattern=path+"*cam2."+htype+".20[23]*.nc"
  print(historyPattern)
  saveFileName=path+caseid+".vars."+htype+"."+year+".nc"
;  saveFileName=path+caseid+".vars."+htype+".2025-39.nc"
  print(saveFileName)
  saveFile=addfile(saveFileName,"c")

;Create UNLIMITED Record Dimension
;  An UNLIMITED record dimension in a netCDF file is critical if the user 
;  plans to ever use the netCDF operators (NCO) to concatenate the file 
;  along the time/record dimension.  
  filedimdef(saveFile,"time",-1,True)

  histVarExtract(historyPattern,saveFile)
end do
