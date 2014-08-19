; These routines create a new vertical structure and then interpolates existing initial
; condition and SAD files to the new strcuture.
;
; NOTE: This routine makes use of ncl and interpic, which must have already been
; built. interpic is a CAM utility.
;
; Author: Mike Mills
; Date: May 23, 2012

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; Generate the initial tracer distribution.


icpath1="/glade/p/work/mmills/inputdata/init/"
icpath2="/glade/p/cesm/wawg/mmills/run/b.e11.B55TRW5CN.f19_g16.MJMtest.013/rest/1989-01-01-00000/"
icpath3="~/inputdata/atm/waccm/ic/"
templatefile="/glade/p/work/mmills/inputdata/init/coords_1.9x2.5_L70.nc"

srcIcWaccm = icpath1+"f_1950-2010_2deg_refc1_wa4_tsmlt.002.cam.i.1987-01-01-00000.nc"
srcIcWACCM5   = icpath2+"b.e11.B55TRW5CN.f19_g16.MJMtest.013.cam.i.1989-01-01-00000.nc"
tmpIcWaccm = icpath3+"tmp.nc"
datestr = systemfunc("date +%y%m%d")
outIcWaccm = icpath3+"b.waccm-tsmlt-mam3_1.9x2.5_L70.cam.i.1987-01-01.c"+datestr+".nc"

print("rm " + outIcWaccm)
system("rm " + outIcWaccm)
print("rm " + tmpIcWaccm)
system("rm " + tmpIcWaccm)

print("interpic -t " + templatefile + " " + srcIcWaccm + " " + outIcWaccm)
system("interpic -t " + templatefile + " " + srcIcWaccm + " " + outIcWaccm)
print(" ")
print(" ")
print("    WACCM file " + srcIcWaccm)
print("    WACCM5 file " + srcIcWACCM5)

; Open the initial conditions file.
WACCMicFile = addfile(outIcWaccm, "w")
WACCM5icFile   = addfile(srcIcWACCM5, "r")

print("rm " + tmpIcWaccm)
system("rm " + tmpIcWaccm)
  
levWACCM = WACCMicFile->lev
levWACCM5   = WACCM5icFile->lev
nlevWACCM5  = dimsizes(levWACCM5)
nlevWACCM  = dimsizes(levWACCM)
print( "nlev WACCM5=" + nlevWACCM5 + ", nlev WACCM=" + nlevWACCM )
if (nlevWACCM5.ne.nlevWACCM) then
  return
end if
  
fields = (/ "OCS", "H2SO4", "S", "SO", "SO3", "SO2", "DMS", "SOAG", "so4_a1", "pom_a1", "soa_a1", "bc_a1", \
 "dst_a1", "ncl_a1", "num_a1", "so4_a2", "soa_a2", "ncl_a2", "num_a2", "dst_a3", \
 "ncl_a3", "so4_a3", "num_a3" /)

nFields = dimsizes(fields)

do i = 0, nFields - 1
  print("copying "+fields(i))
  WACCMicFile->$fields(i)$ = WACCM5icFile->$fields(i)$
end do  
  
print("  IC Done.")

