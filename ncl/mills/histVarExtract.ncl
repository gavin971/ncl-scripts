;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; histVarExtract.ncl
; Used to extract variables from a series of files
;  combining them into 1 file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "ncl-scripts/getHistoryPath.ncl"

undef("histVarExtract")
procedure histVarExtract(historyPattern:string, saveFile:file, varList:string)

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
      
 ; Add default fields to new files, append varList to all
  if .not.(isfilevar(saveFile,"P0")) then
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
  end if
  
  nVars = dimsizes(varList)
  do i = 0, nVars -1
    print("saving "+varList(i)+"...")
    var =addfiles_GetVar(history, historyNames, varList(i))
    saveFile->$varList(i)$ = var
    delete(var)
  end do

end ; procedure histVarExtract


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;   Main
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;caseid="noaur_eco"
;branch="epp"

caseid="125smin"
branch="dust-sulfur"

htype="h0"

path=getHistoryPath(caseid)+"monthly/"

;do year=2027,2040
;  historyPattern=path+"*cam2."+htype+"."+year+"*.nc"
;  historyPattern=path+"*cam2."+htype+".20[23]*.nc"
year=1998
do month=1,4
  if (month.lt.10) then
    monstr="0"+month
  else
    monstr=month
  end if
  historyPattern=path+"all.*cam2."+htype+"."+year+"-"+monstr+"*.nc"
  
  print(historyPattern)
;  saveFileName=path+caseid+".vars."+htype+"."+year+".nc"
  saveFileName=path+caseid+".vars."+htype+"."+year+"-"+monstr+".nc"
  print(saveFileName)
  
; Add variables to existing file, or create new if not existing  
  if (isfilepresent(saveFileName)) then
    saveFile=addfile(saveFileName,"w")
  else
    saveFile=addfile(saveFileName,"c")
;Create UNLIMITED Record Dimension
;  An UNLIMITED record dimension in a netCDF file is critical if the user 
;  plans to ever use the netCDF operators (NCO) to concatenate the file 
;  along the time/record dimension.  
    filedimdef(saveFile,"time",-1,True)
  end if
  
;  varList=(/"NOX","NOY","O3","PS","PSL","T","U","V","OMEGA","Z3","QSUM"/)
  varList=(/"Z3","MIXSULNDBIN","PURSULNDBIN","DUSTNDBIN","SULFRAD",\
            "PURERAD","CORERAD","DUSTRAD"/)

  histVarExtract(historyPattern,saveFile,varList)
  
end do
