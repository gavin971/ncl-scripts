; Sulfur

load "ncl-scripts/core.ncl"
; Generate the initial tracer distribution.

  system("mkdir -p dyn")
  
  srcIcPath = "/home/mills/waccm/inputdata/wa319_smax/waccm319_8_smaxa.cam2.i.2025-01-01-00000.nc"
;  srcIcPath2 = "/home/mills/waccm/inputdata/atm/waccm/ic/carmasadn3.cam2.i.1997-11-01-00000.nc"
  srcIcPath2 = "/home/mills/waccm/inputdata/atm/waccm/ic/nongeo5.cam2.i.1999-01-01-00000.nc"
  icPath    = "dyn/wa3_4x5_dustsulf.pur.extSO2.nogeo.cam2.i.2025-01-01-00000.nc"
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
  
  print("  Adding DUSSUL01 ic... ")
  icFile->DUSSUL01=srcIcFile2->DUSSUL01
  print("  Adding DUSSUL02 ic... ")
  icFile->DUSSUL02=srcIcFile2->DUSSUL02
  print("  Adding DUSSUL03 ic... ")
  icFile->DUSSUL03=srcIcFile2->DUSSUL03
  print("  Adding DUSSUL04 ic... ")
  icFile->DUSSUL04=srcIcFile2->DUSSUL04
  print("  Adding DUSSUL05 ic... ")
  icFile->DUSSUL05=srcIcFile2->DUSSUL05
  print("  Adding DUSSUL06 ic... ")
  icFile->DUSSUL06=srcIcFile2->DUSSUL06
  print("  Adding DUSSUL07 ic... ")
  icFile->DUSSUL07=srcIcFile2->DUSSUL07
  print("  Adding DUSSUL08 ic... ")
  icFile->DUSSUL08=srcIcFile2->DUSSUL08
  print("  Adding DUSSUL09 ic... ")
  icFile->DUSSUL09=srcIcFile2->DUSSUL09
  print("  Adding DUSSUL10 ic... ")
  icFile->DUSSUL10=srcIcFile2->DUSSUL10
  print("  Adding DUSSUL11 ic... ")
  icFile->DUSSUL11=srcIcFile2->DUSSUL11
  print("  Adding DUSSUL12 ic... ")
  icFile->DUSSUL12=srcIcFile2->DUSSUL12
  print("  Adding DUSSUL13 ic... ")
  icFile->DUSSUL13=srcIcFile2->DUSSUL13
  print("  Adding DUSSUL14 ic... ")
  icFile->DUSSUL14=srcIcFile2->DUSSUL14
  print("  Adding DUSSUL15 ic... ")
  icFile->DUSSUL15=srcIcFile2->DUSSUL15
  print("  Adding DUSSUL16 ic... ")
  icFile->DUSSUL16=srcIcFile2->DUSSUL16
  print("  Adding DUSSUL17 ic... ")
  icFile->DUSSUL17=srcIcFile2->DUSSUL17
  print("  Adding DUSSUL18 ic... ")
  icFile->DUSSUL18=srcIcFile2->DUSSUL18
  print("  Adding DUSSUL19 ic... ")
  icFile->DUSSUL19=srcIcFile2->DUSSUL19
  print("  Adding DUSSUL20 ic... ")
  icFile->DUSSUL20=srcIcFile2->DUSSUL20
  print("  Adding DUSSUL21 ic... ")
  icFile->DUSSUL21=srcIcFile2->DUSSUL21
  print("  Adding DUSSUL22 ic... ")
  icFile->DUSSUL22=srcIcFile2->DUSSUL22
  print("  Adding DUSSUL23 ic... ")
  icFile->DUSSUL23=srcIcFile2->DUSSUL23
  print("  Adding DUSSUL24 ic... ")
  icFile->DUSSUL24=srcIcFile2->DUSSUL24
  print("  Adding DUSSUL25 ic... ")
  icFile->DUSSUL25=srcIcFile2->DUSSUL25
  print("  Adding DUSSUL26 ic... ")
  icFile->DUSSUL26=srcIcFile2->DUSSUL26
  print("  Adding DUSSUL27 ic... ")
  icFile->DUSSUL27=srcIcFile2->DUSSUL27
  print("  Adding DUSSUL28 ic... ")
  icFile->DUSSUL28=srcIcFile2->DUSSUL28
  print("  Adding DUSSUL29 ic... ")
  icFile->DUSSUL29=srcIcFile2->DUSSUL29
  print("  Adding DUSSUL20 ic... ")
  icFile->DUSSUL30=srcIcFile2->DUSSUL30
  print("  Adding DUSSUL31 ic... ")
  icFile->DUSSUL31=srcIcFile2->DUSSUL31
  print("  Adding DUSSUL32 ic... ")
  icFile->DUSSUL32=srcIcFile2->DUSSUL32
  print("  Adding DUSSUL33 ic... ")
  icFile->DUSSUL33=srcIcFile2->DUSSUL33
  print("  Adding DUSSUL34 ic... ")
  icFile->DUSSUL34=srcIcFile2->DUSSUL34
  print("  Adding DUSSUL35 ic... ")
  icFile->DUSSUL35=srcIcFile2->DUSSUL35
  print("  Adding DUSSUL36 ic... ")
  icFile->DUSSUL36=srcIcFile2->DUSSUL36
  print("  Adding DUSSUL37 ic... ")
  icFile->DUSSUL37=srcIcFile2->DUSSUL37
  print("  Adding DUSSUL38 ic... ")
  icFile->DUSSUL38=srcIcFile2->DUSSUL38
  
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

