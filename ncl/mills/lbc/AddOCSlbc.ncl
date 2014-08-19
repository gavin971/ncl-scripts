; Sulfur

;load "ncl-scripts/core.ncl"
; Generate the initial tracer distribution.

  datestr = systemfunc("date +%y%m%d")
  srcLbPath = "/glade/p/cesmdata/cseg/inputdata/atm/waccm/lb/"
  srcFile="LBC_1765-2100_1.9x2.5_CCMI_RCP60_za_c130313.nc"
  lbPath    = "/glade/p/work/mmills/inputdata/lb/"
  lbFile="LBC_1765-2100_1.9x2.5_CCMI_RCP60_za_c"+datestr+".nc"
  ; Copy the initial conditions file.
  print("  Creating " + lbPath + lbFile)
  print("    from " + srcLbPath + srcFile)
  system("cp " + srcLbPath + srcFile + " " + lbPath + lbFile)

  ; Open the initial conditions file.
  lbFileH = addfile(lbPath + lbFile, "w")

  ; Add initial conditions for Sulfur
  OCS = lbFileH->CO2_LBC
  OCS@long_name = "OCS lower boundary condition"
  ; OCS: [time ] x [lat | 46]
  OCS( : , :) = 5.12e-10
  print("  Adding OCS lbc... ")
  lbFileH->OCS_LBC = OCS
  print("  LB Done.")

