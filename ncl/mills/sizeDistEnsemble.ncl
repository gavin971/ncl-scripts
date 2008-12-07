; ***********************************************
; sizeDistEnsemble.ncl
; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "tests/historyPath.ncl"
;************************************************

undef("sizeDistEnsemble")
procedure sizeDistEnsemble(caseid)

local

begin
  print("")
  print(caseid)
  historyPath=historyPath(caseid)


end ; procedure sizeDistEnsemble()

caseid="125smin"
sizeDistEnsemble(caseid)
