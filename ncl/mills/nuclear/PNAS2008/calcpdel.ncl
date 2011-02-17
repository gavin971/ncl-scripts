load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

function calcPdel(historyPath:string)

local ps,p0,hyai,hybi,pm

begin

  ; Open the set of history files.
  ;
  ; NOTE: We want to work with the entire set of history files, 
  ; so we can do a time evolution.
  print("  Reading data ...")

  historyNames = systemfunc("ls " + historyPath)
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
  
  pdel=pm(:,1:dimsizes(lev)-1,:,:)-pm(:,0:dimsizes(lev)-2,:,:)
  
  return(pdel)
  
end
