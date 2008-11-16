; Sulfur

load "tests/core.ncl"
; Generate the initial tracer distribution.

  system("mkdir -p dyn")
  
  srcLbPath = "/home/mills/waccm/inputdata/atm/waccm/lb/LBC_Scen=A1b_1950-2050_4x5_sulf.nc"
  lbPath    = "dyn/LBC_Scen=A1b_1950-2050_4x5_OCSx50.nc"
  ; Copy the initial conditions file.
  print("  Creating " + lbPath)
  print("    from " + srcLbPath)
  system("cp " + srcLbPath + " " + lbPath)

  ; Open the initial conditions file.
  lbFile = addfile(lbPath, "w")

  ; Add initial conditions for Sulfur
  OCS = lbFile->OCS_LBC
  ; OCS: [time | 1] x [lev | 66] x [lat | 46] x [lon | 72]
  OCS( : , :, : ) = 2.55e-8
  print("  Adding OCS lbc... ")
  lbFile->OCS_LBC = OCS
  print("  LB Done.")

