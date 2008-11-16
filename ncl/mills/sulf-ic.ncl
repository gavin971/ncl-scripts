; Sulfur Tests
;
; Create some tracer distributions for use with WACCM/CARMA.

WACCM_RES = getenv("WACCM_RES")

if (ismissing(WACCM_RES)) then
;  WACCM_RES  = "2x2.5"
  WACCM_RES  = "4x5"
end if
print("Using " + WACCM_RES + " resolution ...")


load "tests/sulf.ncl"

print("Resolution" + WACCM_RES)

if (WACCM_RES .eq. "10x15") then 
  newic("data/atm/cam/inic/fv/cami_0000-01-01_10x15_L26_c030918.nc")
end if

if (WACCM_RES .eq. "4x5") then 
  newic("data/atm/waccm/ic/wa3_4x5_1950_spinup.cam2.i.1960-01-01-00000.nc")
;  newic("/home/mills/waccm/data/wa319_smax/waccm319_8_smaxa.cam2.i.2025-01-01-00000.nc")
end if

if (WACCM_RES .eq. "2x2.5") then 
  newic("data/atm/waccm/ic/wa3_4x5_1950_spinup.cam2.i.1960-01-01-00000.nc")
end if
