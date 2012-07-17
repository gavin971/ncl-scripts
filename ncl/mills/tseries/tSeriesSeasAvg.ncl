
procedure tSeriesSeasAvg(case:string, season:string, field:string, yearMin:integer, yearMax:integer)

begin
  isWaccm=str_index_of_substr(case,"wcm",0)

  if (ismissing(isWaccm)) then
    rootdir="/data/ccsm4/"
  else
    rootdir="/data/waccm/"
  end if
;  print(rootdir)

  if (season.eq."DJF") then
    months=(/1,2,3/) ; December average given January 1 date, etc.
  else if (season.eq."MAM") then
      months=(/4,5,6/) 
    else if (season.eq."JJA") then
        months=(/7,8,9/) 
      else if (season.eq."SON") then
          months=(/10,11,12/) 
        else
          print("season not supported:"+season)
          exit
        end if
      end if
    end if
  end if

  historyPattern=rootdir+case+"/proc/tseries/monthly/"+case+".cam2.h0zm."+field+"*.nc"
  filename=systemfunc("ls "+historyPattern)
  print("adding file:"+filename(0))

  a=addfile(filename(0),"r")

  print("get time")
  time = a->time
  print("get date")
  date = a->date
  year = date/10000
    
  month = (date-year*10000)/100
  print(" ")
  print("year(0:23) month(0:23)")
  print(year(0:23)+" "+month(0:23))
  
  
  print(" ")
  print("get "+field)
  field = a->$field$
  lat = field&lat
  lev = field&lev
  nlat=dimsizes(lat)
  nlev=dimsizes(lev)
  nyears=yearMax-yearMin+1

  fieldSeasAvg = new((/nyears,nlev,nlat/),"float")
  fieldSeasAvg!0 = "year"
  fieldSeasAvg!1 = "lev"
  fieldSeasAvg!2 = "lat"
  year2=ispan(yearMin,yearMax,1)
  fieldSeasAvg&year = year2
  fieldSeasAvg&lev = lev
  fieldSeasAvg&lat = lat
  
  do y=0,nyears-1
    print("year="+year2(y))
    indmonth = ind(year.eq.year2(y)  .and. month.ge.months(0) .and. month.le.months(2))
    print(indmonth)
    ds=dimsizes(indmonth)
    if (ds(0).ne.3) then
      print("FATAL: ds="+ds)
      exit
    end if
    fieldSeasAvg(y,:,:)=dim_avg(field(lev|:,lat|:,time|indmonth))
  end do
    
  outfile=rootdir+case+"/proc/"+case+".cam2.h0zm."+field+"."+season+"."+yearMin+"-"+yearMax+".nc"
  print("creating "+outfile)
  b=addfile(outfile,"c")
  b->$field$=fieldSeasAvg
  
end

tSeriesSeasAvg("b40.20th.track1.1deg.005", "DJF", "U", 1960, 1979)
