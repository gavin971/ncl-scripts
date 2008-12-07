;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; histVarCombine.ncl
; Used to stitch new variables brought over from Columbia
;  into files previously transfered
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("histVarCombine")
procedure histVarCombine(varsFile:file,vars2File:file)

local varsFile,vars2File

begin
  if (any(varsFile->date.ne.vars2File->date)) then
    print("date misalignment error")
    print(varsFile)
    print(varsFile->date)
    print(vars2File)
    print(vars2File->date)
    exit
  else
    print("dates match")  
  end if
  
  print("combining NOY...")
  varsFile->NOY = vars2File->NOY
  
  print("combining U...")
  varsFile->U = vars2File->U
  
  print("combining V...")
  varsFile->V = vars2File->V
  
  print("combining OMEGA...")
  varsFile->OMEGA = vars2File->OMEGA

end

;Main

caseid="mee01fco"
htype="h3"
path="run/"+caseid+"/"

do year=2028,2049
  varsFileName=path+caseid+".vars."+htype+"."+year+".nc"
  vars2FileName=path+caseid+".vars2."+htype+"."+year+".nc"

  print(varsFileName)
  print(vars2FileName)
  
  varsFile=addfile(varsFileName,"w")
  vars2File=addfile(vars2FileName,"r")

  histVarCombine(varsFile,vars2File)
end do
