; Soot Tests
;
; Create some tracer distributions for use with WACCM/CARMA.

WACCM_RES = getenv("WACCM_RES")
MODEL_PATH = getenv("MODEL_PATH")

if (ismissing(WACCM_RES)) then
  WACCM_RES  = "2x2.5"
;  WACCM_RES  = "4x5"
end if
print("Using " + WACCM_RES + " resolution ...")


load "soot.ncl"

print("Resolution: " + WACCM_RES)

;if (WACCM_RES .eq. "10x15") then 
;  newic("data/atm/cam/inic/fv", "cami_0000-01-01_10x15_L26_c030918.nc", "output.nc")
;end if

;if (WACCM_RES .eq. "4x5") then 
;;  newic("data/atm/waccm/ic/wa3_4x5_1950_spinup.cam2.i.1960-01-01-00000.nc")
;  newic("/home/mills/waccm/data/wa319_smax", \
;        "waccm319_8_smaxa.cam2.i.2025-01-01-00000.nc", \
;        "waccm319_8_smaxa_soot5Tg.cam2.i.2025-01-01-00000.nc", \
;        1.8031e-5)
;end if

if (WACCM_RES .eq. "2x2.5") then 
;  newic("data/atm/waccm/ic/wa3_4x5_1950_spinup.cam2.i.1960-01-01-00000.nc")

;  newic(MODEL_PATH+"/waccm/inputdata/atm/waccm/ic", \
;       "f2000.e10r02.2deg.waccm.005.cam2.i.0017-01-01-00000.nc", \
;       "f40.2000.2deg.waccm.carma.bc5_0tg.cam2.i.0017-01-01-00000.nc", \
;       4.50775e-06)

;  newic(MODEL_PATH+"/waccm/inputdata/atm/waccm/ic", \
;       "b40.2000.2deg.wcm.carma.bc.002.0Tg.cam2.i.0001-05-15-00000.nc", \
;       "b40.2000.2deg.wcm.carma.bc.002.5Tg_4x5.cam2.i.0001-05-15-00000.nc", \
;       5.0, 4.0, 5.0)

;  newic(MODEL_PATH+"/waccm/inputdata/atm/waccm/ic", \
;       "b40.2000.2deg.wcm.carma.bc.002.0Tg.cam2.i.0001-05-15-00000.nc", \
;       "b40.2000.2deg.wcm.carma.bc.002.5Tg_8x10.cam2.i.0001-05-15-00000.nc", \
;       5.0, 8.0, 10.0)

  newic(MODEL_PATH+"/waccm/inputdata/atm/waccm/ic", \
       "b40.2000.2deg.wcm.carma.bc.002.0Tg.cam2.i.0001-05-15-00000.nc", \
       "b40.2000.2deg.wcm.carma.bc.002.5Tg.cam2.i.0001-05-15-00000.nc", \
       5.0)

end if
