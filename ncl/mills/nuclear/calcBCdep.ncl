load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;        float a2x_BCPHIDRY(time, lat, lon) ;
;                a2x_BCPHIDRY:units = "kg/m2/s" ;
;                a2x_BCPHIDRY:long_name = "drydep of hydrophilic black carbon" ;
;                a2x_BCPHIDRY:cell_methods = "time: mean" ;
;        float a2x_BCPHIWET(time, lat, lon) ;
;                a2x_BCPHIWET:units = "kg/m2/s" ;
;                a2x_BCPHIWET:long_name = "wetdep of hydrophilic black carbon" ;
;                a2x_BCPHIWET:cell_methods = "time: mean" ;
;        float a2x_BCPHODRY(time, lat, lon) ;
;                a2x_BCPHODRY:units = "kg/m2/s" ;
;                a2x_BCPHODRY:long_name = "drydep of hydrophobic black carbon" ;
;                a2x_BCPHODRY:cell_methods = "time: mean" ;
;        float a2x_OCPHIDRY(time, lat, lon) ;
;                a2x_OCPHIDRY:units = "kg/m2/s" ;
;                a2x_OCPHIDRY:long_name = "drydep of hydrophilic organic carbon" ;
;                a2x_OCPHIDRY:cell_methods = "time: mean" ;
;        float a2x_OCPHIWET(time, lat, lon) ;
;                a2x_OCPHIWET:units = "kg/m2/s" ;
;                a2x_OCPHIWET:long_name = "wetdep of hydrophilic organic carbon" ;
;                a2x_OCPHIWET:cell_methods = "time: mean" ;
;        float a2x_OCPHODRY(time, lat, lon) ;
;                a2x_OCPHODRY:units = "kg/m2/s" ;
;                a2x_OCPHODRY:long_name = "drydep of hydrophobic organic carbon" ;
;                a2x_OCPHODRY:cell_methods = "time: mean" ;

  ; Some constants we will need:
  pi      = 4.0 * atan(1.0)
  dg2rad  = pi / 180.0   ; rad deg^-1
  g       = 9.806               ; m s^-2 
  re      = 6.37122e06          ; m

;testpath1="/data5/mmills/b40.rcp4_5.2deg.wcm.carma.bc5tgIndPak.002/h3"
;cntlpath1="/data5/mmills/b40.rcp4_5.2deg.wcm.carma.bc0tg.001/h3"
testPath2="/data5/mmills/b40.rcp4_5.2deg.wcm.carma.bc5tgIP.002/h3/"
cntlPath2="/data5/mmills/b40.rcp4_5.2deg.wcm.carma.bc0tg.002/h3/"

  testNames = testPath2+systemfunc("ls " + testPath2)
  print("Reading files ...")
  print((/ testNames /))
  tests = addfiles(testNames, "r")
  test1 = addfile(testNames(0),"r")
  ListSetType(tests, "cat")
  
  print("reading test a2x_BCPHIDRY")
  bcphidry_t2 = addfiles_GetVar(tests, testNames, "a2x_BCPHIDRY")
  print("reading test a2x_BCPHIWET")
  bcphiwet_t2 = addfiles_GetVar(tests, testNames, "a2x_BCPHIWET")
  print("reading test a2x_BCPHODRY")
  bcphodry_t2 = addfiles_GetVar(tests, testNames, "a2x_BCPHODRY")
  print("reading test SOOT01DD")
  soot01dry_t2 = addfiles_GetVar(tests, testNames, "SOOT01DD")
  print("reading test SOOT01SFWET")
  soot01wet_t2 = addfiles_GetVar(tests, testNames, "SOOT01SFWET")
;  ocphidry_t2 = addfiles_GetVar(tests, testNames, "a2x_OCPHIDRY")
;  ocphiwet_t2 = addfiles_GetVar(tests, testNames, "a2x_OCPHIWET")
;  ocphodry_t2 = addfiles_GetVar(tests, testNames, "a2x_OCPHODRY")
  date_t2     = addfiles_GetVar(tests, testNames, "date")
  time     = addfiles_GetVar(tests, testNames, "time")
  gw          = test1->gw
  lat         = test1->lat
  lon         = test1->lon
  
  cntlNames = cntlPath2+systemfunc("ls " + cntlPath2)
;  cntlNames = str_sub_str(testNames,".bc5tgIP.",".bc0tg.")
  print("Reading files ...")
  print((/ cntlNames /))
  cntls = addfiles(cntlNames, "r")
  ListSetType(cntls, "cat")
  
  print("reading cntl a2x_BCPHIDRY")
  bcphidry_c2 = addfiles_GetVar(cntls, cntlNames, "a2x_BCPHIDRY")
  print("reading cntl a2x_BCPHIWET")
  bcphiwet_c2 = addfiles_GetVar(cntls, cntlNames, "a2x_BCPHIWET")
  print("reading cntl a2x_BCPHODRY")
  bcphodry_c2 = addfiles_GetVar(cntls, cntlNames, "a2x_BCPHODRY")
;  ocphidry_c2 = addfiles_GetVar(cntls, cntlNames, "a2x_OCPHIDRY")
;  ocphiwet_c2 = addfiles_GetVar(cntls, cntlNames, "a2x_OCPHIWET")
;  ocphodry_c2 = addfiles_GetVar(cntls, cntlNames, "a2x_OCPHODRY")
  date_c2     = addfiles_GetVar(cntls, cntlNames, "date")
  
  nDates = dimsizes(date_t2)
  d1t2=ind(date_t2.eq."20120515")
  d1c2=ind(date_c2.eq."20120515")
  d2c2=ind(date_c2.eq.date_t2(nDates-1))
  print("20120515: "+d1t2+" "+d1c2+" "+d2c2)
  
  ; Determine some of the dimensions.
  nLats = dimsizes(lat) 
  nLons = dimsizes(lon)
    
  wgt = new((/ nLats, nLons /), double)
  wgt = conform(wgt, gw*2.*pi/nLons*re*re*24.*3600./1.e9, 0) ; m2*(s/day)*(Tg/kg)
  
  bcphidry_t2day = new(nDates,double)
  bcphidry_t2day!0="time"
  bcphidry_t2day&time=time
  bcphidry_t2day@units="Tg/day"
  
  bcphiwet_t2day=bcphidry_t2day
  bcphodry_t2day=bcphidry_t2day  
  soot01wet_t2day=bcphidry_t2day
  soot01dry_t2day=bcphidry_t2day
  
  bcphidry_c2day=bcphidry_t2day
  bcphiwet_c2day=bcphidry_t2day
  bcphodry_c2day=bcphidry_t2day
  
  bcphidry_t2day(:) = wgt_areasum2(bcphidry_t2, wgt, 0) ; Tg/day
  bcphiwet_t2day(:) = wgt_areasum2(bcphiwet_t2, wgt, 0) ; Tg/day
  bcphodry_t2day(:) = wgt_areasum2(bcphodry_t2, wgt, 0) ; Tg/day
  soot01dry_t2day(:) = wgt_areasum2(soot01dry_t2, wgt, 0) ; Tg/day
  soot01wet_t2day(:) = wgt_areasum2(soot01wet_t2, wgt, 0) ; Tg/day
  
  bcphidry_c2day(:) = wgt_areasum2(bcphidry_c2(d1c2:d2c2,:,:), wgt, 0) ; Tg/day
  bcphiwet_c2day(:) = wgt_areasum2(bcphiwet_c2(d1c2:d2c2,:,:), wgt, 0) ; Tg/day
  bcphodry_c2day(:) = wgt_areasum2(bcphodry_c2(d1c2:d2c2,:,:), wgt, 0) ; Tg/day
  print(bcphidry_t2day+"	"+bcphiwet_t2day+"	"+bcphodry_t2day+\
  "	"+bcphidry_c2day+"	"+bcphiwet_c2day+"	"+bcphodry_c2day)
;  print(bcphidry_t2day)
  
  bcphidry_t2day@long_name="drydep of hydrophylic black carbon, 5Tg global"
  bcphiwet_t2day@long_name="wetdep of hydrophylic black carbon, 5Tg global"
  bcphodry_t2day@long_name="drydep of hydrophobic black carbon, 5Tg global"
  soot01dry_t2day@long_name="drydep of CARMA black carbon, 5Tg global"
  soot01wet_t2day@long_name="wetdep of CARMA black carbon, 5Tg global"
  bcphidry_c2day@long_name="drydep of hydrophylic black carbon, control global"
  bcphiwet_c2day@long_name="wetdep of hydrophylic black carbon, control global"
  bcphodry_c2day@long_name="drydep of hydrophobic black carbon, control global"
  
  outPath="/data5/mmills/b40.rcp4_5.2deg.wcm.carma.bc5tgIP.002/proc/"
  outFile="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.002.bcdep.nc"
  test5 = systemfunc("mkdir -p " + outPath)
  test6 = systemfunc("rm " + outPath+outFile)
  outfile = addfile(outPath+outFile,"c")
  outfile->date=date_t2
  outfile->bcphidry_test=bcphidry_t2day
  outfile->bcphiwet_test=bcphiwet_t2day
  outfile->bcphodry_test=bcphodry_t2day
  outfile->soot01dry_test=soot01dry_t2day
  outfile->soot01wet_test=soot01wet_t2day
  
  outfile->bcphidry_cntl=bcphidry_c2day
  outfile->bcphiwet_cntl=bcphiwet_c2day
  outfile->bcphodry_cntl=bcphodry_c2day
  

