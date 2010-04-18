load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

historyRun="vgsoot08"
controlRun="vgcontrol02"
;rootPath="/data5/mmills"
;rootPath="/Volumes/Data/Models/waccm/carma-soot/run"
rootPath="/Data/Models/waccm/run/carma-soot"

  print((/historyRun/))
;  saveMeanVars(historyRun,True)
;  saveMeanLOxFiles(historyRun)
  historyPath = rootPath+"/"+historyRun+"/h0zm/" 
  historyPattern = historyPath+historyRun+".cam2.h0zm.[01][0-9].nc"
  historyNames = systemfunc("ls " + historyPattern)
;  print((/historyNames/))
  historySDPattern = historyPath+historyRun+".cam2.h0zmsd.[01][0-9].nc"
  historySDNames = systemfunc("ls " + historySDPattern)
    
  a=addfile(historyNames(0),"r")
  lat=a->lat
  b=addfile(historySDNames(0),"r")
  O3_col_Jan=a->O3_col
  O3_colSD_Jan=b->O3_col(0,:)
  TS_Jan=a->TS
  TSSD_Jan=b->TS(0,:)
  
;  print("  Reading data ...")
  JanPattern=historyPath+historyRun+".cam2.h0zm.*-01.nc"
  print(JanPattern)
  historyNamesJan = systemfunc("ls " + JanPattern)
  history = addfiles(historyNamesJan, "r")
  ListSetType(history, "cat")
  
  ; Read in the fields we will need.
  O3_col_Jans = addfiles_GetVar(history, historyNames, "O3_col")
  print(dimsizes(O3_col_Jans))
  print(O3_col_Jans!0)
  print(O3_col_Jans!1)
  
  O3_col_Jans_mean = dim_avg(O3_col_Jans(lat|:,case|:))
  print(dimsizes(O3_col_Jans_mean))
  O3_col_Jans_var = dim_variance(O3_col_Jans(lat|:,case|:))
  print(dimsizes(O3_col_Jans_mean))
 
  fmt = "%8.3f"
  fmt2 = "%5.0f"
  print(sprintf(fmt2, lat)+" "+sprintf(fmt, O3_col_Jan)+" "+sprintf(fmt, O3_colSD_Jan)+\
     " "+sprintf(fmt, O3_col_Jans_mean)+" "+sprintf(fmt, O3_col_Jans_var)+" "+\
     sprintf(fmt, O3_colSD_Jan^2))
  
  TS_Jans = addfiles_GetVar(history, historyNames, "TS")
  TS_Jans_mean = dim_avg(TS_Jans(lat|:,case|:))
  TS_Jans_var = dim_variance(TS_Jans(lat|:,case|:))
  
  print(sprintf(fmt2, lat)+" "+sprintf(fmt, TS_Jan)+" "+sprintf(fmt, TSSD_Jan)+\
     " "+sprintf(fmt, TS_Jans_mean)+" "+sprintf(fmt, TS_Jans_var)+" "+\
     sprintf(fmt, TSSD_Jan^2))
  
  
