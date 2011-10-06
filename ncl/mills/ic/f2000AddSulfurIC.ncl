; Sulfur

;load "../core.ncl"
; Generate the initial tracer distribution.

  modelpath=getenv("MODEL_PATH")
  icpath=modelpath+"/ccsm/inputdata/atm/waccm/ic/"
  
  srcIcPath  = icpath+"f40.2000.track1.4deg.001.cam2.i.0013-01-01-00000.nc"
  srcIcPath2 = icpath+"nogeo400SO2.cam2.i.1999-01-01-00000.nc"
  icPath     = icpath+"f40.2000.track1.4deg.001.sulfgas.cam2.i.0013-01-01-00000.nc"
  
  ; Copy the initial conditions file.
  print("  Creating " + icPath)
  print("    from " + srcIcPath)
  system("cp " + srcIcPath + " " + icPath)

  ; Open the initial conditions file.
  icFile = addfile(icPath, "w")
  srcIcFile2 = addfile(srcIcPath2, "r")

  ; Add initial conditions for Sulfur
  print("  Adding OCS ic... ")
  icFile->OCS=srcIcFile2->OCS
  print("  Adding S ic... ")
  icFile->S=srcIcFile2->S
  print("  Adding SO ic... ")
  icFile->SO=srcIcFile2->SO
  print("  Adding SO2 ic... ")
  icFile->SO2=srcIcFile2->SO2
  print("  Adding SO3 ic... ")
  icFile->SO3=srcIcFile2->SO3
  print("  Adding HSO3 ic... ")
  icFile->HSO3=srcIcFile2->HSO3
  print("  Adding H2SO4 ic... ")
  icFile->H2SO4=srcIcFile2->H2SO4

  print("  IC Done.")

