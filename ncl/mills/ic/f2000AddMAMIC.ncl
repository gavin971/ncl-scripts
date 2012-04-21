load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; Generate the initial tracer distribution.

;modelpath=getenv("MODEL_PATH")
;icpath=modelpath+"/ccsm/inputdata/atm/waccm/ic/"
;icpath2=modelpath+"/ccsm/inputdata/atm/cam/inic/fv/"
;srcIcPath  = icpath+"f2000.e10r02.2deg.waccm.005.cam2.i.0017-01-01-00000.nc"
;srcIcPath2 = icpath2+"cami-mam3_0000-01-01_1.9x2.5_L30_c090306.nc"
;icPath     = icpath+"f2000.waccm-mam3_1.9x2.5_cam2.i.0017-01-01.c120215.nc"

icpath="/glade/proj3/cseg/inputdata/atm/waccm/ic/"
icpath2="/glade/proj3/cseg/inputdata/atm/cam/inic/fv/"
srcIcPath  = icpath+"f2000.waccm-mam3_1.9x2.5_L70.cam2.i.0017-01-01.c120410.nc"
srcIcPath2 = icpath2+"trop_strat_mam3_chem_2000-01-01_1.9x2.5_L30_c120222.nc"
icPath     = icpath+"f2000.waccm-mam3_1.9x2.5_L70_SOA.cam2.i.0017-01-01.c120420.nc"

; Copy the initial conditions file.
print("  Creating " + icPath)
print("    from " + srcIcPath)
system("cp " + srcIcPath + " " + icPath)

; Open the initial conditions file.
WACCMicFile = addfile(icPath, "w")
CAMicFile   = addfile(srcIcPath2, "r")
  
levWACCM = WACCMicFile->lev
levCAM   = CAMicFile->lev
nlevCAM  = dimsizes(levCAM)
nlevWACCM  = dimsizes(levWACCM)
print( "nlev CAM=" + nlevCAM + ", nlev WACCM=" + nlevWACCM )
  
fields = (/ "H2SO4", "SO2", "DMS", "SOAG", "SOA", "so4_a1", "pom_a1", "soa_a1", "bc_a1", \
 "dst_a1", "ncl_a1", "num_a1", "so4_a2", "soa_a2", "ncl_a2", "num_a2", "dst_a3", \
 "ncl_a3", "so4_a3", "num_a3" /)

nFields = dimsizes(fields)

do i = 0, nFields - 1
  print("interpolating "+fields(i))
  specCAM = CAMicFile->$fields(i)$(time|:,lat|:,lon|:,lev|:)
  specWACCM = linint1_Wrap(levCAM, specCAM(time|:,lat|:,lon|:,lev|:), False, levWACCM, 0)
  print ("max=" + max(specWACCM) + ", min=" + min(specWACCM))
; set values at levels above CAM top to a small number
  specWACCM( time|:, lat|:, lon|:, {lev|0:levCAM(0)} ) = 1.0d-45
; set values at lowest level to those at CAM lowest level  
  specWACCM( time|:, lat|:, lon|:, lev|nlevWACCM - 1 ) = specCAM(time|:, lat|:, lon|:, lev|nlevCAM - 1)
  WACCMicFile->$fields(i)$ = specWACCM(time|:,lev|:,lat|:,lon|:)
end do  
  
print("  IC Done.")

