load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

rootpath="/glade/p/work/mmills/inputdata/emis/"
datestr = systemfunc("date +%y%m%d")

ox = getenv("ox")
if (.not.ismissing(ox)) then
  if ((ox.eq."OCS").or.(ox.eq."CH4")) then
    print("ox = "+ox)
    ox_chml=ox+"_CHML"

    if (ox.eq."OCS") then
      file1850     ="SO2_emission_OCS_oxidation_elev_1850cyclical_WACCM5_c150224.nc"
      file1850_1857="SO2_emission_OCS_oxidation_elev_1850-1857_WACCM5_c150227.nc"
      file1883_2000="SO2_emission_OCS_oxidation_elev_1883-2000_WACCM5_c150228.nc"
      file1985_2010="SO2_emission_OCS_oxidation_elev_1985-2010_WACCM5_c150301.nc"
      outfile      ="SO2_emission_OCS_oxidation_elev_1849-2101_WACCM5_c"+datestr+".nc"
      factor=1.0
      long_name="SO2 production derived from OCS chemical loss rate"
    end if

    if (ox.eq."CH4") then
      file1850     ="CH4_oxidation_2H2O_emission_elev_1850_WACCM5_c150228.nc"
      file1850_1857="CH4_oxidation_2H2O_emission_elev_1850-1857_WACCM5_c150227.nc"
      file1883_2000="CH4_oxidation_2H2O_emission_elev_1883-2000_WACCM5_c150228.nc"
      file1985_2010="CH4_oxidation_2H2O_emission_elev_1985-2010_WACCM5_c150301.nc"
      outfile      ="H2O_emission_CH4_oxidationx2_elev_1849-2101_WACCM5_c"+datestr+".nc"
      factor=2.0
      long_name="H2O production derived from 2x CH4 chemical loss rate"
    end if

    print("opening "+file1850)
    ncid1850=addfile(rootpath+file1850,"r")
    print("opening "+file1850_1857)
    ncid1850_1857=addfile(rootpath+file1850_1857,"r")
    print("opening "+file1883_2000)
    ncid1883_2000=addfile(rootpath+file1883_2000,"r")
    print("opening "+file1985_2010)
    ncid1985_2010=addfile(rootpath+file1985_2010,"r")

    lev=ncid1850->lev
    lat=ncid1850->lat
    lon=ncid1850->lon
    nlev=dimsizes(lev)
    nlat=dimsizes(lat)
    nlon=dimsizes(lon)

    cmd="rm "+rootpath+outfile
    print((/cmd/))
    a=systemfunc(cmd)
    print("creating "+rootpath+outfile)
    ncidout=addfile(rootpath+outfile,"c")

    ; Create an UNLIMITED record dimension in the output netCDF file.  This is critical if 
    ;  the user plans to ever use NCO to concatenate the file along the time/record dimension.
    filedimdef(ncidout,"time",-1,True)
    ;
    ; For a nicer looking netCDF, create a "new line" character.
    ; This is not necessary.
    ;
    nl = integertochar(10)  ; newline character
    ;
    ; Define global attributes.
    ;
    ; globalAtt can be of any type. Here logical is used by convention.
    ;
    globalAtt             = True

    globalAtt@creation_date     = nl+\
        systemfunc("date")

    globalAtt@creator = nl+\
        "Mike Mills, mmills@ucar.edu"
    globalAtt@title = nl+\
        long_name ;
    globalAtt@Sources = nl+\
         "Derived from WACCM runs covering the following time periods:"+nl+\
         "  1849-1850: b.e13.B1850W5TCN.f19_g16.igwsProgAer.noMEGAN.008.rhminl.8900, monthly averages years 1-32"+nl+\
         "  1851-1857: b.e13.B20TRW5TCN.f19_g16.beta13_stratmodal.001"+nl+\
         "  1858-1882: interpolated corresponding months from 1857 to 1883"+nl+\
         "  1883-1984: b.e13.B20TRW5TCN.f19_g16.beta13_stratmodal.002"+nl+\
         "  1985-2010: b.e13.B20TRW5TCN.f19_g16.beta13_stratmodal.006.0-15N"+nl+\
         "  2011-2101: repeat year 2010"
    globalAtt@history = nl+\
         "Interpolated and concatenated by program ExtFrcCH4-H2O_OCS-SO2oxidation.ncl"

    fileattdef( ncidout, globalAtt )

    ncidout->lev=lev
    ncidout->lat=lat
    ncidout->lon=lon

    ncidout->hyam=ncid1850->hyam
    ncidout->hybm=ncid1850->hybm
    ncidout->hyai=ncid1850->hyai
    ncidout->hybi=ncid1850->hybi
    ncidout->P0=ncid1850->P0
    ncidout->ilev=ncid1850->ilev
    ncidout->gw=ncid1850->gw

    CHML1850=ncid1850->$ox_chml$
    CHML1850_1857=ncid1850_1857->$ox_chml$
    CHML1883_2000=ncid1883_2000->$ox_chml$
    CHML1985_2010=ncid1985_2010->$ox_chml$


    date=yyyymm_time(1849, 2101, "integer")
    date=date*100+15 ; day in mid-month
    date(1::12)=date(1::12)-1 ; 14th February
    date@long_name = "current date (YYYYMMDD)"
    delete(date&time)
;    time=yyyymmdd_to_yyyyfrac(date,0)
;    time=floattoint((time-1849)*365)
;    time@units="days since January 1, 1849"
    
    print("writing date")
    ncidout->date=date

    ntime=dimsizes(date)
    print("creating CHML")
    CHML=new((/ ntime, nlev, nlat, nlon /),"float")
    CHML@long_name=long_name
    CHML@units="/cm3/s"
    CHML!0="time"
    CHML!1="lev"
    CHML!2="lat"
    CHML!3="lon"
;    CHML&time=date&time
    CHML&lev=lev
    CHML&lat=lat
    CHML&lon=lon

    print("1850...")
    CHML(0:11,:,:,:)=(/CHML1850/)
    CHML(12:23,:,:,:)=(/CHML1850/)

    print("1850-1857...")
    m11=(1851-1849)*12
    m12=(1858-1849)*12-1
    m21=(1851-1850)*12
    m22=(1858-1850)*12-1
    CHML(m11:m12,:,:,:)=(/CHML1850_1857(m21:m22,:,:,:)/)

    print("1883-1984...")
    m11=(1883-1849)*12
    m12=(1985-1849)*12-1
    m21=(1883-1883)*12
    m22=(1985-1883)*12-1
    CHML(m11:m12,:,:,:)=(/CHML1883_2000(m21:m22,:,:,:)/)

    print("1985-2010...")
    m11=(1985-1849)*12
    m12=(2011-1849)*12-1
    m21=(1985-1985)*12
    m22=(2011-1985)*12-1
    CHML(m11:m12,:,:,:)=(/CHML1985_2010(m21:m22,:,:,:)/)

    m2010=(2010-1849)*12
    print("2011-2101 (copying 2010)...")
    do year = 2011,2101
      m1=(year-1849)*12
      CHML(m1:m1+11,:,:,:)=(/CHML(m2010:m2010+11,:,:,:)/)
    end do

    m1857=(1857-1849)*12
    m1883=(1883-1849)*12
    ny=1883.-1857.

    print("1858-1882 (interpolating)...")
    do year=1858,1882
      m1=(year-1849)*12
      f1=(year-1857.)/ny
      f2=(1883.-year)/ny
      CHML(m1:m1+11,:,:,:)=f2*CHML(m1857:m1857+11,:,:,:)+f1*CHML(m1883:m1883+11,:,:,:)
    end do

    if (factor.ne.1.0) then
      print("multiplying by "+factor)
      CHML=factor*CHML
    end if

    print("writing "+ox_chml)
    ncidout->$ox_chml$=CHML
    print("closing "+outfile)
    delete(ncidout)
  else
    print("Set environment variable ox to either OCS or CH4, then try again")
  end if
else
  print("Set environment variable ox to either OCS or CH4, then try again")
end if
