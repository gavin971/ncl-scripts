
procedure tSeriesSeasAvg(case:string, season:string, field:string, yearMin:integer, yearMax:integer)

begin
  host=systemfunc("hostname")
  print(host)

  if (host.eq."cmip5") then
    isWaccm=str_index_of_substr(case,"wcm",0)
    if (ismissing(isWaccm)) then
       rootdir="/data/ccsm4/"
    else
       rootdir="/data/waccm/"
    end if
  else
    rootdir="/Volumes/MillsPassport/cesm/case/"
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
  fieldarr = a->$field$
  lat = fieldarr&lat
  lev = fieldarr&lev
  nlat=dimsizes(lat)
  nlev=dimsizes(lev)
  nyears=yearMax-yearMin+1

  fieldSeasAvg = new((/nyears,nlev,nlat/),"float")
  fieldSeasAvg!0 = "time"
  fieldSeasAvg!1 = "lev"
  fieldSeasAvg!2 = "lat"
  year2=ispan(yearMin,yearMax,1)
  year2!0="time"
  date2=year2*10000+months(0)*100+1
  date2!0="time"
;  fieldSeasAvg&year = year2
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
    fieldSeasAvg(y,:,:)=dim_avg(fieldarr(lev|:,lat|:,time|indmonth))
  end do
    
  outfile=rootdir+case+"/proc/"+case+".cam2.h0zm."+field+"."+season+"."+yearMin+"-"+yearMax+".nc"
  if (isfilepresent(outfile)) then
    print("deleting existing file "+outfile)
    delfile=systemfunc("rm "+outfile)
    print(delfile)
  end if
  print("creating "+outfile)
  b=addfile(outfile,"c")
  b->year=year2
  b->date=date2
  b->P0=a->P0
  b->gw=a->gw
  b->hyai=a->hyai
  b->hyam=a->hyam
  b->hybi=a->hybi
  b->hybm=a->hybm
  b->slat=a->slat
  b->slon=a->slon  
  b->$field$=fieldSeasAvg
  
end

do i = 5,9
  tSeriesSeasAvg("b40.20th.track1.1deg.00"+i, "DJF", "U", 1960, 1979)
  tSeriesSeasAvg("b40.20th.track1.1deg.00"+i, "DJF", "U", 1986, 2005)
end do
