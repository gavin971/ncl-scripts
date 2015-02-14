load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

  
  ostype=getenv("OSTYPE")
  host=getenv("HOSTNAME")
  
  rootpath="/glade/p/work/mmills/data/SAGE_II_v7.00/CCMI_Filled_Data_V2/"
  if (ostype.eq."darwin") then
    rootpath="/Volumes/Data/Data/SAGE2/SAGE_II_v7.00/CCMI_Filled_Data_V2/"
  end if
  
  sagedir=rootpath+"ASCII/"
  infiles=systemfunc("ls -1 "+sagedir)
  ntimes=dimsizes(infiles)
  print("ntimes="+ntimes)

  lat=(/-77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5, -42.5, -37.5, -32.5, -27.5, -22.5, \
        -17.5, -12.5,  -7.5,  -2.5, 2.5, 7.5,  12.5,  17.5,  22.5,  27.5,  32.5,  37.5, \
         42.5,  47.5,  52.5,  57.5,  62.5,  67.5,  72.5,  77.5 /)
  lat@long_name = "latitude"
  lat@units = "degrees_north"

  nlat=dimsizes(lat)
  ncol = 41
  nalt=70
  nhead=9
  ntot=nhead+nalt

  ext1020=new((/ntimes,nalt,nlat/),"float")
  ext1020!0="time"
  ext1020!1="alt"
  ext1020!2="lat"
  ext1020&lat=lat
  ext1020@_FillValue=-1.0
  ext1020@long_name="Extinction at 1020 nm"
  ext1020@units="km-1"
  ext525=ext1020
  ext525@long_name="Extinction at 525 nm"
  ext452=ext1020
  ext452@long_name="Extinction at 452 nm"
  ext386=ext1020
  ext386@long_name="Extinction at 386 nm"
  
  date=new(ntimes,"float")
  
  do i=0,ntimes-1
  
    filename=infiles(i)
    print(i+" of "+ntimes+" "+filename)
    
    filechar=stringtochar(filename)
    yearmonthstr=chartostring(filechar(0:5))
;    print("yearmonthstr="+yearmonthstr)
    yearmonth=stringtoint(yearmonthstr)
    date(i)=yyyymm_to_yyyyfrac(yearmonth,0.5)
    print("  "+date(i))

    numRows = numAsciiRow(sagedir+filename)
;    print("numRows="+numRows)

    do l=0,nlat-1
      print("      "+l+" "+lat(l))

      data = readAsciiTable(sagedir+filename, ncol, "float", (/nhead+l*ntot,numRows-(l+1)*ntot/))
      if (l.eq.0) then
        alt=data(:,0)
        alt@long_name="Altitude"
        alt@units="km"
;        print(alt)
        ext1020&alt=alt
        ext525&alt=alt
        ext452&alt=alt
        ext386&alt=alt     
      end if

      ext1020(i,:,l)=(/data(:,1)/)
      ext525(i,:,l)=(/data(:,6)/)
      ext452(i,:,l)=(/data(:,11)/)
      ext386(i,:,l)=(/data(:,16)/)

    end do
    
  end do
  ext1020&time=date
  ext525&time=date
  ext452&time=date
  ext386&time=date
  
  a=systemfunc("mkdir "+rootpath+"netCDF")
  b=addfile(rootpath+"netCDF/CCMI_Filled_Data_V2.nc","c")
  b->ext1020=ext1020
  b->ext525=ext525
  b->ext452=ext452
  b->ext386=ext386
  delete(b)
  
  
