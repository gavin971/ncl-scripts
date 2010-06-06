;runname="b40.20th.track1.2deg.wcm.007"
;wadir="/data3/mmills/"+runname+"/h0zm.annual"
;year1=1885
;year2=1912

print("in fixdate.ncl:")
runname=getenv("runname")
print(runname)
wadir=getenv("anndir")
print(wadir)
year1=stringtoint(getenv("year1"))
print(year1)
year2=stringtoint(getenv("year2"))
print(year2)

do year = year1, year2
  date=year*10000+701
  if ( year .lt. 10 ) then
      yearstr="000"+year
  else
    if ( year .lt. 100 ) then
      yearstr="00"+year
    else
      if ( year .lt. 1000 ) then
        yearstr="0"+year
      else
        yearstr=""+year
      end if
    end if
  end if
  
  filename=wadir+"/"+runname+".cam2.h0zm."+yearstr+".annual.nc"
  print(date+"  "+filename)
  a=addfile(filename,"w")
  a->date=(/date/)
  delete(a)
end do
  
