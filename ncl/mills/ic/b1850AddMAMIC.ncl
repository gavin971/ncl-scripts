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


ostype=getenv("OSTYPE")
host=getenv("HOSTNAME")

if (ostype.eq."darwin") then
  icpath1="/Volumes/MillsPassport/cesm/ic/"
  icpath2=icpath1
  icpath3=icpath1
else
  if(str_index_of_substr(host,"mirage",-1).eq.0)
    icpath1="/glade/proj3/cseg/inputdata/ccsm4_init/b40.1850.track1.2deg.wcm.007/0156-01-01/"
    icpath2="/glade/proj3/cseg/inputdata/ccsm4_init/b40_1850_2d_r07c5cn_160jp/0070-01-01/"
    icpath3="/glade/proj3/cseg/inputdata/atm/waccm/ic/"
    templatefile="/glade/proj2/acd/fvitt/waccm-ic-L70/coords_1.9x2.5_L70.nc"
  else if (host.eq."modeling1") then
    icpath1="/data2/mmills/"  
    icpath2=icpath1
    icpath3=icpath1
  else
    icpath1="/data/waccm/"  
    icpath2=icpath1
    icpath3=icpath1
  end if
end if
srcIcWaccm = icpath1+"b40.1850.track1.2deg.wcm.007.cam2.i.0156-01-01-00000.nc"
srcIcCam   = icpath2+"b40_1850_2d_r07c5cn_160jp.cam.i.0070-01-01-00000.nc"
outIcWaccm = icpath3+"b1850.waccm-mam3_1.9x2.5_L70.cam2.i.0156-01-01.c120523.nc"

print("  Creating " + outIcWaccm)
print("    from WACCM4 file " + srcIcWaccm)
system("interpic -t " + templatefile + " " + srcIcWaccm + " " + outIcWaccm)

print("    and CAM5 file " + srcIcCam)

; Open the initial conditions file.
WACCMicFile = addfile(outIcWaccm, "w")
CAMicFile   = addfile(srcIcCam, "r")
  
levWACCM = WACCMicFile->lev
levCAM   = CAMicFile->lev
nlevCAM  = dimsizes(levCAM)
nlevWACCM  = dimsizes(levWACCM)
print( "nlev CAM=" + nlevCAM + ", nlev WACCM=" + nlevWACCM )
  
fields = (/ "H2SO4", "SO2", "DMS", "SOAG", "so4_a1", "pom_a1", "soa_a1", "bc_a1", \
 "dst_a1", "ncl_a1", "num_a1", "so4_a2", "soa_a2", "ncl_a2", "num_a2", "dst_a3", \
 "ncl_a3", "so4_a3", "num_a3" /)

nFields = dimsizes(fields)

do i = 0, nFields - 1
  print("interpolating "+fields(i))
  specCAM = CAMicFile->$fields(i)$(time|:,lat|:,lon|:,lev|:)
  specWACCM = linint1_Wrap(levCAM, specCAM(time|:,lat|:,lon|:,lev|:), False, levWACCM, 0)
  print ("max=" + max(specWACCM) + ", min=" + min(specWACCM))
; set values at levels above CAM top to a small number
  specWACCM( time|:, lat|:, lon|:, {lev|0:levCAM(0)} ) = 1.0d-36
; set values at lowest level to those at CAM lowest level  
  specWACCM( time|:, lat|:, lon|:, lev|nlevWACCM - 1 ) = specCAM(time|:, lat|:, lon|:, lev|nlevCAM - 1)
  WACCMicFile->$fields(i)$ = specWACCM(time|:,lev|:,lat|:,lon|:)
end do  
  
print("  IC Done.")

