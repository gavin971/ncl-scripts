; Add OCS to WACCM 1850 LBC

rootpathic="/glade/p/work/mmills/inputdata/init/"
WACCM_1850IC_Filename="b.e12.B1850W5CN.f19_g16.002.cam.i.0200-01-01-00000.nc"

WACCM_1850IC_File=addfile(rootpathic+WACCM_1850IC_Filename,"r")

creation_date = systemfunc("date +%y%m%d")
outFilename="b.e12.B1850W5CN.f19_g16.002.cam.i.0200-01-01_OCS_c"+creation_date+".nc"
cmd="cp "+rootpathic+WACCM_1850IC_Filename+" "+rootpathic+outFilename
print(cmd)
cp=systemfunc(cmd)

outFile=addfile(rootpathic+outFilename,"w")

icpath2="/glade/p/cesm/wawg/mmills/run/b.e11.B55TRW5CN.f19_g16.MJMtest.013/rest/1989-01-01-00000/"
OCS_ICfilename   = "b.e11.B55TRW5CN.f19_g16.MJMtest.013.cam.i.1989-01-01-00000.nc"
OCS_ICfile=addfile(icpath2+OCS_ICfilename,"r")

OCSratio=337.5/512.0
print("Multiplying OCS by "+OCSratio+"...")
OCS=OCS_ICfile->OCS
OCS=OCS*OCSratio
OCS@history="1989 values scaled to 1850 by factor of 337.5/512.0="+OCSratio
outFile->OCS=OCS

;
; For a nicer looking netCDF, create a "new line" character.
; This is not necessary.
;
  nl = integertochar(10)  ; newline character
  
outFile@history=nl+\
          "Created from "+WACCM_1850IC_Filename+nl+\
          "using program b1850addSulfurIC.ncl"+nl+\
          "OCS from "+OCS_ICfilename+nl+\
          "OCS scaled to 1850 by factor of 337.5/512.0="+OCSratio
outFile@creation_date     = nl+\               
        systemfunc("date")
outFile@creator = nl+\
      "Mike Mills, mmills@ucar.edu"
      
print("Done")
delete(outFile)
