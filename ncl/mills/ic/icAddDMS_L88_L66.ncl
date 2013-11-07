; create new IC files for 88 and 66 vertical levels by adding DMS from  
; a CAM4 run
;
; NOTE: This routine makes use of ncl and interpic, which must have already been
; built. interpic is a CAM utility.
;
; Author: Mike Mills
; Date: Nov 6, 2013

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


file88_in  = "/glade/p/cesmdata/cseg/inputdata/atm/waccm/ic/" + \
   "f40.2008-2011.2deg.wcm.sd.carma.sulf.002.cam2.i.2010-01-01-00000.nc"


fileDMS_in = "/glade/p/acd/ccmi/f_1975_2010_refc1sd_CESM1_CAM4chem.001/rest/2000-01-01-00000/" + \
   "f_1975_2010_refc1sd_CESM1_CAM4chem.001.cam.i.2000-01-01-00000.nc"

icPath = "/glade/p/work/mmills/inputdata/init/"
file88_out=icPath+"f19.wcm.sd.carma.sulf.DMS.L88.cam2.i.2010-01-01-00000.nc"
file66_out=icPath+"f19.wcm.carma.sulf.DMS.L66.cam2.i.2010-01-01-00000.nc"
filetemp  =icPath+"temp.nc"

print ("rm " + filetemp)
system("rm " + filetemp)

print ("interpic -t " + file88_in + " " + fileDMS_in + " " + filetemp)
system("interpic -t " + file88_in + " " + fileDMS_in + " " + filetemp)
print(" ")
print(" ")
print ("cp " + file88_in + " " + file88_out)
system("cp " + file88_in + " " + file88_out)
print(" ")
print(" ")

a=addfile(file88_out,"w")
b=addfile(filetemp,"r")

print("copying DMS from "+filetemp+" to "+file88_out)
a->DMS=b->DMS

template66="/glade/p/cesmdata/cseg/inputdata/atm/waccm/ic/cami_2000-07-01_1.9x2.5_L66_c040928.nc"
print ("rm " + file66_out)
system("rm " + file66_out)
print ("interpic -t " + template66 + " " + file88_out + " " + file66_out)
system("interpic -t " + template66 + " " + file88_out + " " + file66_out)

