; Extend WACCM LBC file back to 850AD

rootpath="/glade/p/work/mmills/inputdata/lb/"
LM_ann_Filename="ghg_pmip3_850-2007_annual_c100517.v2.nc"
; vars: CO2, N2O, CH4, f11, f12
; date: 8500701 to 20070701 annual

rootpathic="/glade/p/work/mmills/inputdata/init/"
WACCM_1850IC_Filename="b.e12.B1850W5CN.f19_g16.002.cam.i.0200-01-01-00000.nc"

LM_ann_File=addfile(rootpath+LM_ann_Filename,"r")

ind850=ind(LM_ann_File->date.eq.8500701)
CO2_850=(LM_ann_File->CO2(ind850))
CH4_850=(LM_ann_File->CH4(ind850))
N2O_850=(LM_ann_File->N2O(ind850))
print("ind850:"+ind850+" CO2_850="+CO2_850+" CH4_850="+CH4_850+" N2O_850="+N2O_850)

ind1850=ind(LM_ann_File->date.eq.18500701)
CO2_1850=(LM_ann_File->CO2(ind1850))
CH4_1850=(LM_ann_File->CH4(ind1850))
N2O_1850=(LM_ann_File->N2O(ind1850))
print("ind1850:"+ind1850+" CO2_1850="+CO2_1850+" CH4_1850="+CH4_1850+" N2O_1850="+N2O_1850)

WACCM_1850IC_File=addfile(rootpathic+WACCM_1850IC_Filename,"r")

creation_date = systemfunc("date +%y%m%d")
outFilename="b.e12.B1850W5CN.f19_g16.002.cam.i.0200-01-01_850ADscaled_c"+creation_date+".nc"
cmd="cp "+rootpathic+WACCM_1850IC_Filename+" "+rootpathic+outFilename
print(cmd)
cp=systemfunc(cmd)

outFile=addfile(rootpathic+outFilename,"w")

CO2ratio=CO2_850/CO2_1850
print("Multiplying CO2 by "+CO2ratio+"...")
CO2=outFile->CO2
CO2=CO2*CO2ratio
CO2@history="1850 values scaled to 850AD by factor of "+CO2ratio
outFile->CO2=CO2

CH4ratio=CH4_850/CH4_1850
print("Multiplying CH4 by "+CH4ratio+"...")
CH4=outFile->CH4
CH4=CH4*CH4ratio
CH4@history="1850 values scaled to 850AD by factor of "+CH4ratio
outFile->CH4=CH4

N2Oratio=N2O_850/N2O_1850
print("Multiplying N2O by "+N2Oratio+"...")
N2O=outFile->N2O
N2O=N2O*N2Oratio
N2O@history="1850 values scaled to 850AD by factor of "+N2Oratio
outFile->N2O=N2O

;
; For a nicer looking netCDF, create a "new line" character.
; This is not necessary.
;
  nl = integertochar(10)  ; newline character
  
outFile@history=nl+\
          "Created from "+WACCM_1850IC_Filename+nl+\
          "using program LastMilleniumIC.ncl"+nl+\
          "1850 values for CO2, CH4, and N2O scaled to 850AD values"+nl+\
          "CO2 scaled by factor of "+CO2ratio+nl+\
          "CH4 scaled by factor of "+CH4ratio+nl+\
          "N2O scaled by factor of "+N2Oratio
outFile@creation_date     = nl+\               
        systemfunc("date")
outFile@creator = nl+\
      "Mike Mills, mmills@ucar.edu"
      
print("Done")
delete(outFile)
