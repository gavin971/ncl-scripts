;runname="b40.20th.track1.2deg.wcm.007"
;wadir="/data3/mmills/"+runname+"/h0zm.annual"
;year1=1885
;year2=1912

print("in fixdate-W5.ncl:")
runname="b.e11.B1850W5CN.f19_g16.007"
print("runname="+(/runname/))
wadir="/glade/proj3/CESM/wawg/runs/b.e11.B1850W5CN.f19_g16.007a/h0"
print("wadir="+(/wadir/))
year1=1
print("year1="+(/year1/))
year_end=21
print("year_end="+(/year_end/))
searchstr="cam.h0"
print("searchstr="+(/searchstr/))

do year = year1, year_end
  do month = 1,12
    if (month .eq. 12) then
      month2=1
      year2=100+year+1
    else
      month2=month+1
;      year2=100+year
      year2=year
    end if
    
    date=year2*10000+month2*100+01
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
    
    if ( month .lt. 10 ) then
      monthstr="0"+month
    else
      monthstr=""+month
    end if
    

    filename=wadir+"/"+runname+"."+searchstr+"."+yearstr+"-"+monthstr+".nc"
    print(date+"  "+filename)
    a=addfile(filename,"w")
    a->date=(/date/)
;    a->time=a->time+36500
;    a->time=a->time-36500
    delete(a)
  end do
end do
  
