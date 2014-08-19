; Sulfur

load "ncl-scripts/core.ncl"
; Generate the initial tracer distribution.

  system("mkdir -p dyn")
  
  srcIcPath = "/home/mills/waccm/inputdata/atm/waccm/ic/wa3_4x5_dustsulf.extpur2.cam2.i.1996-01-01-00000.nc"
  icPath    = "dyn/wa3_4x5_dustsulf.extpur2.ocsx50.cam2.i.1996-01-01-00000.nc"
  ; Copy the initial conditions file.
  print("  Creating " + icPath)
  print("    from " + srcIcPath)
  system("cp " + srcIcPath + " " + icPath)

  ; Open the initial conditions file.
  icFile = addfile(icPath, "w")

  ; Add initial conditions for Sulfur
  OCS = icFile->OCS
  lev = icFile->lev
  ; OCS: [time | 1] x [lev | 66] x [lat | 46] x [lon | 72]
  OCS( 0, 55: , :, : ) = OCS( 0, 55: , :, : )*50.
  OCS( 0, 54 , :, : )  = OCS( 0, 54 , :, : )*40.
  OCS( 0, 53 , :, : ) = OCS( 0, 53 , :, : )*30.
  OCS( 0, 52 , :, : ) = OCS( 0, 52 , :, : )*20.
  OCS( 0, 51 , :, : ) = OCS( 0, 51 , :, : )*10.
  print("  Adding OCS ic... ")
  icFile->OCS = OCS
  print("  IC Done.")

