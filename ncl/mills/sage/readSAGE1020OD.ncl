load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("readSAGE1020OD")
function readSAGE1020OD ()

begin
  ostype=getenv("OSTYPE")
  host=getenv("HOSTNAME")
  
  sagedir="/glade/p/work/mmills/data/SAGE_II_v7.00/"
  if (ostype.eq."darwin") then
    sagedir="/Volumes/Data/Data/SAGE2/SAGE_II_v7.00/"
  end if
    
  sagefile="SAGE_1020_OD_Filled.dat"

  lat=(/-77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5, -42.5, -37.5, -32.5, -27.5, -22.5, \
        -17.5, -12.5,  -7.5,  -2.5, 2.5, 7.5,  12.5,  17.5,  22.5,  27.5,  32.5,  37.5, \
         42.5,  47.5,  52.5,  57.5,  62.5,  67.5,  72.5,  77.5 /)
  lat@long_name = "latitude"
  lat@units = "degrees_north"

  nlat=dimsizes(lat)
  ncol=nlat+1

  data = readAsciiTable(sagedir+sagefile, ncol, "float", 2)

  date=data(:,0)
  date@long_name = "year"

  ndate=dimsizes(date)

  od=(/data(:,1:)/)
  od@_FillValue = 9.999
  od=10.^od
  od@_FillValue = -999
  od!0="time"
  od&time=date
  od!1="lat"
  od&lat=lat
  
  printVarSummary(od)

  return(od)
end

