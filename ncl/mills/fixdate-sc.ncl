;runname="b40.20th.track1.2deg.wcm.007"
;wadir="/data3/mmills/"+runname+"/h0zm.annual"
;year1=1885
;year2=1912

print("in fixdate-sc.ncl:")
forcingfile=getenv("forcingfile")
print("forcingfile="+(/forcingfile/))
year1=stringtoint(getenv("year1"))
print("year1="+(/year1/))
year2=stringtoint(getenv("year2"))
print("year2="+(/year2/))

a=addfile(forcingfile,"w")
date=a->date
i=0

do year = year1, year2
  do month = 1, 12
    if (month .eq. 2) then
      date(i)=year*10000+month*100+14
    else
      date(i)=year*10000+month*100+15
    end if
    i=i+1
  end do 
end do

; Adjust first and last date to cover full range of model run:
date(0) = year1*10000+101
date(i-1) = (year2+1)*10000+102

a->date=(/date/)
delete(a)
  
