; Sulfur

load "tests/core.ncl"
; Generate the initial tracer distribution.
undef("newic")
procedure newic(srcIcPath:string)

begin
  system("mkdir -p dyn")

  icPath    = "dyn/wa3_4x5_1950_spinup_sulfur.cam2.i.1960-01-01-00000.nc"
  ; Copy the initial conditions file.
  print("  Creating " + icPath)
  print("    from " + srcIcPath)
  system("cp " + srcIcPath + " " + icPath)

  ; Open the initial conditions file.
  icFile = addfile(icPath, "w")

  ; Add initial conditions for Sulfur
  tmp = icFile->CFC11
  ; OCS
  tmp( 0, : , :, : ) = 510.e-12 * 60.08 / 28.964 ; OCS mmr = (510 pptv) * M_OCS / M_air
  tmp@long_name="OCS"
  print("  Adding OCS ... ")
  icFile->OCS = tmp
  ; Set rest of the sulfur species to small number
  tmp( 0, : , :, : ) = 1.e-25 
  tmp@long_name="S"
  print("  Adding S ... ")
  icFile->S = tmp
  tmp@long_name="SO"
  print("  Adding SO ... ")
  icFile->SO = tmp
  tmp@long_name="SO2_g"
  print("  Adding SO2 ... ")
  icFile->SO2_g = tmp
  tmp@long_name="SO3"
  print("  Adding SO3 ... ")
  icFile->SO3 = tmp
  tmp@long_name="HSO3"
  print("  Adding HSO3 ... ")
  icFile->HSO3 = tmp
  tmp@long_name="H2SO4"
  print("  Adding H2SO4 ... ")
  icFile->H2SO4 = tmp
  print("  Done.")
end

