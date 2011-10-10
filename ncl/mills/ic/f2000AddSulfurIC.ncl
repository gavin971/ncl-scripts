; Sulfur

;load "../core.ncl"
; Generate the initial tracer distribution.

  modelpath=getenv("MODEL_PATH")
  icpath=modelpath+"/ccsm/inputdata/atm/waccm/ic/"
  
  srcIcPath  = icpath+"f40.2000.track1.4deg.001.cam2.i.0013-01-01-00000.nc"
  srcIcPath2 = icpath+"nogeo400SO2.cam2.i.1999-01-01-00000.nc"
  icPath     = icpath+"f40.2000.track1.4deg.001.sulfbins.cam2.i.0013-01-01-00000.nc"
  
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

  print("  Adding PURSUL01 ic... ")
  icFile->PURSUL01=srcIcFile2->PURSUL01
  print("  Adding PURSUL02 ic... ")
  icFile->PURSUL02=srcIcFile2->PURSUL02
  print("  Adding PURSUL03 ic... ")
  icFile->PURSUL03=srcIcFile2->PURSUL03
  print("  Adding PURSUL04 ic... ")
  icFile->PURSUL04=srcIcFile2->PURSUL04
  print("  Adding PURSUL05 ic... ")
  icFile->PURSUL05=srcIcFile2->PURSUL05
  print("  Adding PURSUL06 ic... ")
  icFile->PURSUL06=srcIcFile2->PURSUL06
  print("  Adding PURSUL07 ic... ")
  icFile->PURSUL07=srcIcFile2->PURSUL07
  print("  Adding PURSUL08 ic... ")
  icFile->PURSUL08=srcIcFile2->PURSUL08
  print("  Adding PURSUL09 ic... ")
  icFile->PURSUL09=srcIcFile2->PURSUL09
  print("  Adding PURSUL10 ic... ")
  icFile->PURSUL10=srcIcFile2->PURSUL10
  print("  Adding PURSUL11 ic... ")
  icFile->PURSUL11=srcIcFile2->PURSUL11
  print("  Adding PURSUL12 ic... ")
  icFile->PURSUL12=srcIcFile2->PURSUL12
  print("  Adding PURSUL13 ic... ")
  icFile->PURSUL13=srcIcFile2->PURSUL13
  print("  Adding PURSUL14 ic... ")
  icFile->PURSUL14=srcIcFile2->PURSUL14
  print("  Adding PURSUL15 ic... ")
  icFile->PURSUL15=srcIcFile2->PURSUL15
  print("  Adding PURSUL16 ic... ")
  icFile->PURSUL16=srcIcFile2->PURSUL16
  print("  Adding PURSUL17 ic... ")
  icFile->PURSUL17=srcIcFile2->PURSUL17
  print("  Adding PURSUL18 ic... ")
  icFile->PURSUL18=srcIcFile2->PURSUL18
  print("  Adding PURSUL19 ic... ")
  icFile->PURSUL19=srcIcFile2->PURSUL19
  print("  Adding PURSUL20 ic... ")
  icFile->PURSUL20=srcIcFile2->PURSUL20
  print("  Adding PURSUL21 ic... ")
  icFile->PURSUL21=srcIcFile2->PURSUL21
  print("  Adding PURSUL22 ic... ")
  icFile->PURSUL22=srcIcFile2->PURSUL22
  print("  Adding PURSUL23 ic... ")
  icFile->PURSUL23=srcIcFile2->PURSUL23
  print("  Adding PURSUL24 ic... ")
  icFile->PURSUL24=srcIcFile2->PURSUL24
  print("  Adding PURSUL25 ic... ")
  icFile->PURSUL25=srcIcFile2->PURSUL25
  print("  Adding PURSUL26 ic... ")
  icFile->PURSUL26=srcIcFile2->PURSUL26
  print("  Adding PURSUL27 ic... ")
  icFile->PURSUL27=srcIcFile2->PURSUL27
  print("  Adding PURSUL28 ic... ")
  icFile->PURSUL28=srcIcFile2->PURSUL28
  print("  Adding PURSUL29 ic... ")
  icFile->PURSUL29=srcIcFile2->PURSUL29
  print("  Adding PURSUL20 ic... ")
  icFile->PURSUL30=srcIcFile2->PURSUL30
  print("  Adding PURSUL31 ic... ")
  icFile->PURSUL31=srcIcFile2->PURSUL31
  print("  Adding PURSUL32 ic... ")
  icFile->PURSUL32=srcIcFile2->PURSUL32
  print("  Adding PURSUL33 ic... ")
  icFile->PURSUL33=srcIcFile2->PURSUL33
  print("  Adding PURSUL34 ic... ")
  icFile->PURSUL34=srcIcFile2->PURSUL34
  print("  Adding PURSUL35 ic... ")
  icFile->PURSUL35=srcIcFile2->PURSUL35
  print("  Adding PURSUL36 ic... ")
  icFile->PURSUL36=srcIcFile2->PURSUL36
  print("  Adding PURSUL37 ic... ")
  icFile->PURSUL37=srcIcFile2->PURSUL37
  print("  Adding PURSUL38 ic... ")
  icFile->PURSUL38=srcIcFile2->PURSUL38

  print("  IC Done.")

