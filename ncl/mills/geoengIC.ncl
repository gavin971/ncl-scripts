; Sulfur

load "tests/core.ncl"
; Generate the initial tracer distribution.

  system("mkdir -p dyn")
  
  srcIcPath = "/home/mills/waccm/inputdata/wa319_smax/waccm319_8_smaxa.cam2.i.2025-01-01-00000.nc"
  srcIcPath2 = "/home/mills/waccm/inputdata/atm/waccm/ic/geoengn4.cam2.i.1998-03-01-00000.nc"
  icPath    = "dyn/wa3_4x5_dustsulf.nopur.extSO2.ocsx50.cam2.i.2025-01-01-00000.nc"
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
  
  print("  Adding DUSTSUL01 ic... ")
  icFile->DUSTSUL01=srcIcFile2->DUSTSUL01
  print("  Adding DUSTSUL02 ic... ")
  icFile->DUSTSUL02=srcIcFile2->DUSTSUL02
  print("  Adding DUSTSUL03 ic... ")
  icFile->DUSTSUL03=srcIcFile2->DUSTSUL03
  print("  Adding DUSTSUL04 ic... ")
  icFile->DUSTSUL04=srcIcFile2->DUSTSUL04
  print("  Adding DUSTSUL05 ic... ")
  icFile->DUSTSUL05=srcIcFile2->DUSTSUL05
  print("  Adding DUSTSUL06 ic... ")
  icFile->DUSTSUL06=srcIcFile2->DUSTSUL06
  print("  Adding DUSTSUL07 ic... ")
  icFile->DUSTSUL07=srcIcFile2->DUSTSUL07
  print("  Adding DUSTSUL08 ic... ")
  icFile->DUSTSUL08=srcIcFile2->DUSTSUL08
  print("  Adding DUSTSUL09 ic... ")
  icFile->DUSTSUL09=srcIcFile2->DUSTSUL09
  print("  Adding DUSTSUL10 ic... ")
  icFile->DUSTSUL10=srcIcFile2->DUSTSUL10
  print("  Adding DUSTSUL11 ic... ")
  icFile->DUSTSUL11=srcIcFile2->DUSTSUL11
  print("  Adding DUSTSUL12 ic... ")
  icFile->DUSTSUL12=srcIcFile2->DUSTSUL12
  print("  Adding DUSTSUL13 ic... ")
  icFile->DUSTSUL13=srcIcFile2->DUSTSUL13
  print("  Adding DUSTSUL14 ic... ")
  icFile->DUSTSUL14=srcIcFile2->DUSTSUL14
  print("  Adding DUSTSUL15 ic... ")
  icFile->DUSTSUL15=srcIcFile2->DUSTSUL15
  print("  Adding DUSTSUL16 ic... ")
  icFile->DUSTSUL16=srcIcFile2->DUSTSUL16
  print("  Adding DUSTSUL17 ic... ")
  icFile->DUSTSUL17=srcIcFile2->DUSTSUL17
  print("  Adding DUSTSUL18 ic... ")
  icFile->DUSTSUL18=srcIcFile2->DUSTSUL18
  print("  Adding DUSTSUL19 ic... ")
  icFile->DUSTSUL19=srcIcFile2->DUSTSUL19
  print("  Adding DUSTSUL20 ic... ")
  icFile->DUSTSUL20=srcIcFile2->DUSTSUL20
  print("  Adding DUSTSUL21 ic... ")
  icFile->DUSTSUL21=srcIcFile2->DUSTSUL21
  print("  Adding DUSTSUL22 ic... ")
  icFile->DUSTSUL22=srcIcFile2->DUSTSUL22
  print("  Adding DUSTSUL23 ic... ")
  icFile->DUSTSUL23=srcIcFile2->DUSTSUL23
  print("  Adding DUSTSUL24 ic... ")
  icFile->DUSTSUL24=srcIcFile2->DUSTSUL24
  print("  Adding DUSTSUL25 ic... ")
  icFile->DUSTSUL25=srcIcFile2->DUSTSUL25
  print("  Adding DUSTSUL26 ic... ")
  icFile->DUSTSUL26=srcIcFile2->DUSTSUL26
  print("  Adding DUSTSUL27 ic... ")
  icFile->DUSTSUL27=srcIcFile2->DUSTSUL27
  print("  Adding DUSTSUL28 ic... ")
  icFile->DUSTSUL28=srcIcFile2->DUSTSUL28
  print("  Adding DUSTSUL29 ic... ")
  icFile->DUSTSUL29=srcIcFile2->DUSTSUL29
  print("  Adding DUSTSUL20 ic... ")
  icFile->DUSTSUL30=srcIcFile2->DUSTSUL30
  print("  Adding DUSTSUL31 ic... ")
  icFile->DUSTSUL31=srcIcFile2->DUSTSUL31
  print("  Adding DUSTSUL32 ic... ")
  icFile->DUSTSUL32=srcIcFile2->DUSTSUL32
  print("  Adding DUSTSUL33 ic... ")
  icFile->DUSTSUL33=srcIcFile2->DUSTSUL33
  print("  Adding DUSTSUL34 ic... ")
  icFile->DUSTSUL34=srcIcFile2->DUSTSUL34
  print("  Adding DUSTSUL35 ic... ")
  icFile->DUSTSUL35=srcIcFile2->DUSTSUL35
  print("  Adding DUSTSUL36 ic... ")
  icFile->DUSTSUL36=srcIcFile2->DUSTSUL36
  print("  Adding DUSTSUL37 ic... ")
  icFile->DUSTSUL37=srcIcFile2->DUSTSUL37
  print("  Adding DUSTSUL38 ic... ")
  icFile->DUSTSUL38=srcIcFile2->DUSTSUL38
  
  print("  Adding DUST01 ic... ")
  icFile->DUST01=srcIcFile2->DUST01
  print("  Adding DUST02 ic... ")
  icFile->DUST02=srcIcFile2->DUST02
  print("  Adding DUST03 ic... ")
  icFile->DUST03=srcIcFile2->DUST03
  print("  Adding DUST04 ic... ")
  icFile->DUST04=srcIcFile2->DUST04
  print("  Adding DUST05 ic... ")
  icFile->DUST05=srcIcFile2->DUST05
  print("  Adding DUST06 ic... ")
  icFile->DUST06=srcIcFile2->DUST06
  print("  Adding DUST07 ic... ")
  icFile->DUST07=srcIcFile2->DUST07
  print("  Adding DUST08 ic... ")
  icFile->DUST08=srcIcFile2->DUST08
  print("  Adding DUST09 ic... ")
  icFile->DUST09=srcIcFile2->DUST09
  print("  Adding DUST10 ic... ")
  icFile->DUST10=srcIcFile2->DUST10
  print("  Adding DUST11 ic... ")
  icFile->DUST11=srcIcFile2->DUST11
  print("  Adding DUST12 ic... ")
  icFile->DUST12=srcIcFile2->DUST12
  print("  Adding DUST13 ic... ")
  icFile->DUST13=srcIcFile2->DUST13
  print("  Adding DUST14 ic... ")
  icFile->DUST14=srcIcFile2->DUST14
  print("  Adding DUST15 ic... ")
  icFile->DUST15=srcIcFile2->DUST15
  print("  Adding DUST16 ic... ")
  icFile->DUST16=srcIcFile2->DUST16
  print("  Adding DUST17 ic... ")
  icFile->DUST17=srcIcFile2->DUST17
  print("  Adding DUST18 ic... ")
  icFile->DUST18=srcIcFile2->DUST18
  print("  Adding DUST19 ic... ")
  icFile->DUST19=srcIcFile2->DUST19
  print("  Adding DUST20 ic... ")
  icFile->DUST20=srcIcFile2->DUST20
  print("  Adding DUST21 ic... ")
  icFile->DUST21=srcIcFile2->DUST21
  print("  Adding DUST22 ic... ")
  icFile->DUST22=srcIcFile2->DUST22
  print("  Adding DUST23 ic... ")
  icFile->DUST23=srcIcFile2->DUST23
  print("  Adding DUST24 ic... ")
  icFile->DUST24=srcIcFile2->DUST24
  print("  Adding DUST25 ic... ")
  icFile->DUST25=srcIcFile2->DUST25
  print("  Adding DUST26 ic... ")
  icFile->DUST26=srcIcFile2->DUST26
  print("  Adding DUST27 ic... ")
  icFile->DUST27=srcIcFile2->DUST27
  print("  Adding DUST28 ic... ")
  icFile->DUST28=srcIcFile2->DUST28
  print("  Adding DUST29 ic... ")
  icFile->DUST29=srcIcFile2->DUST29
  print("  Adding DUST20 ic... ")
  icFile->DUST30=srcIcFile2->DUST30
  print("  Adding DUST31 ic... ")
  icFile->DUST31=srcIcFile2->DUST31
  print("  Adding DUST32 ic... ")
  icFile->DUST32=srcIcFile2->DUST32
  print("  Adding DUST33 ic... ")
  icFile->DUST33=srcIcFile2->DUST33
  print("  Adding DUST34 ic... ")
  icFile->DUST34=srcIcFile2->DUST34
  print("  Adding DUST35 ic... ")
  icFile->DUST35=srcIcFile2->DUST35
  print("  Adding DUST36 ic... ")
  icFile->DUST36=srcIcFile2->DUST36
  print("  Adding DUST37 ic... ")
  icFile->DUST37=srcIcFile2->DUST37
  print("  Adding DUST38 ic... ")
  icFile->DUST38=srcIcFile2->DUST38

  print("  Adding MIXSUL01 ic... ")
  icFile->MIXSUL01=srcIcFile2->MIXSUL01
  print("  Adding MIXSUL02 ic... ")
  icFile->MIXSUL02=srcIcFile2->MIXSUL02
  print("  Adding MIXSUL03 ic... ")
  icFile->MIXSUL03=srcIcFile2->MIXSUL03
  print("  Adding MIXSUL04 ic... ")
  icFile->MIXSUL04=srcIcFile2->MIXSUL04
  print("  Adding MIXSUL05 ic... ")
  icFile->MIXSUL05=srcIcFile2->MIXSUL05
  print("  Adding MIXSUL06 ic... ")
  icFile->MIXSUL06=srcIcFile2->MIXSUL06
  print("  Adding MIXSUL07 ic... ")
  icFile->MIXSUL07=srcIcFile2->MIXSUL07
  print("  Adding MIXSUL08 ic... ")
  icFile->MIXSUL08=srcIcFile2->MIXSUL08
  print("  Adding MIXSUL09 ic... ")
  icFile->MIXSUL09=srcIcFile2->MIXSUL09
  print("  Adding MIXSUL10 ic... ")
  icFile->MIXSUL10=srcIcFile2->MIXSUL10
  print("  Adding MIXSUL11 ic... ")
  icFile->MIXSUL11=srcIcFile2->MIXSUL11
  print("  Adding MIXSUL12 ic... ")
  icFile->MIXSUL12=srcIcFile2->MIXSUL12
  print("  Adding MIXSUL13 ic... ")
  icFile->MIXSUL13=srcIcFile2->MIXSUL13
  print("  Adding MIXSUL14 ic... ")
  icFile->MIXSUL14=srcIcFile2->MIXSUL14
  print("  Adding MIXSUL15 ic... ")
  icFile->MIXSUL15=srcIcFile2->MIXSUL15
  print("  Adding MIXSUL16 ic... ")
  icFile->MIXSUL16=srcIcFile2->MIXSUL16
  print("  Adding MIXSUL17 ic... ")
  icFile->MIXSUL17=srcIcFile2->MIXSUL17
  print("  Adding MIXSUL18 ic... ")
  icFile->MIXSUL18=srcIcFile2->MIXSUL18
  print("  Adding MIXSUL19 ic... ")
  icFile->MIXSUL19=srcIcFile2->MIXSUL19
  print("  Adding MIXSUL20 ic... ")
  icFile->MIXSUL20=srcIcFile2->MIXSUL20
  print("  Adding MIXSUL21 ic... ")
  icFile->MIXSUL21=srcIcFile2->MIXSUL21
  print("  Adding MIXSUL22 ic... ")
  icFile->MIXSUL22=srcIcFile2->MIXSUL22
  print("  Adding MIXSUL23 ic... ")
  icFile->MIXSUL23=srcIcFile2->MIXSUL23
  print("  Adding MIXSUL24 ic... ")
  icFile->MIXSUL24=srcIcFile2->MIXSUL24
  print("  Adding MIXSUL25 ic... ")
  icFile->MIXSUL25=srcIcFile2->MIXSUL25
  print("  Adding MIXSUL26 ic... ")
  icFile->MIXSUL26=srcIcFile2->MIXSUL26
  print("  Adding MIXSUL27 ic... ")
  icFile->MIXSUL27=srcIcFile2->MIXSUL27
  print("  Adding MIXSUL28 ic... ")
  icFile->MIXSUL28=srcIcFile2->MIXSUL28
  print("  Adding MIXSUL29 ic... ")
  icFile->MIXSUL29=srcIcFile2->MIXSUL29
  print("  Adding MIXSUL20 ic... ")
  icFile->MIXSUL30=srcIcFile2->MIXSUL30
  print("  Adding MIXSUL31 ic... ")
  icFile->MIXSUL31=srcIcFile2->MIXSUL31
  print("  Adding MIXSUL32 ic... ")
  icFile->MIXSUL32=srcIcFile2->MIXSUL32
  print("  Adding MIXSUL33 ic... ")
  icFile->MIXSUL33=srcIcFile2->MIXSUL33
  print("  Adding MIXSUL34 ic... ")
  icFile->MIXSUL34=srcIcFile2->MIXSUL34
  print("  Adding MIXSUL35 ic... ")
  icFile->MIXSUL35=srcIcFile2->MIXSUL35
  print("  Adding MIXSUL36 ic... ")
  icFile->MIXSUL36=srcIcFile2->MIXSUL36
  print("  Adding MIXSUL37 ic... ")
  icFile->MIXSUL37=srcIcFile2->MIXSUL37
  print("  Adding MIXSUL38 ic... ")
  icFile->MIXSUL38=srcIcFile2->MIXSUL38

  print("  IC Done.")

