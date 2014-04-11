load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

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
datesec=a->datesec
time_bnds=a->time_bnds
CO2=a->CO2
H=a->H
NO=a->NO
O=a->O
O2=a->O2
O3=a->O3
PS=a->PS

date_2=array_append_record(date(0:0),date,0)
datesec_2=array_append_record(datesec(0:0),datesec,0)
time_bnds_2=array_append_record(time_bnds(0:0,:),time_bnds,0)
CO2_2=array_append_record(CO2(0:0,:,:),CO2,0)
H_2=array_append_record(H(0:0,:,:),H,0)
NO_2=array_append_record(NO(0:0,:,:),NO,0)
O_2=array_append_record(O(0:0,:,:),O,0)
O2_2=array_append_record(O2(0:0,:,:),O2,0)
O3_2=array_append_record(O3(0:0,:,:),O3,0)
PS_2=array_append_record(PS(0:0,:),PS,0)

nTimes=dimsizes(date_2)

date_3=array_append_record(date_2,date_2(nTimes-1:nTimes-1),0)
datesec_3=array_append_record(datesec_2,datesec_2(nTimes-1:nTimes-1),0)
time_bnds_3=array_append_record(time_bnds_2,time_bnds_2(nTimes-1:nTimes-1,:),0)
CO2_3=array_append_record(CO2_2,CO2_2(nTimes-1:nTimes-1,:,:),0)
H_3=array_append_record(H_2,H_2(nTimes-1:nTimes-1,:,:),0)
NO_3=array_append_record(NO_2,NO_2(nTimes-1:nTimes-1,:,:),0)
O_3=array_append_record(O_2,O_2(nTimes-1:nTimes-1,:,:),0)
O2_3=array_append_record(O2_2,O2_2(nTimes-1:nTimes-1,:,:),0)
O3_3=array_append_record(O3_2,O3_2(nTimes-1:nTimes-1,:,:),0)
PS_3=array_append_record(PS_2,PS_2(nTimes-1:nTimes-1,:),0)


date_3(0)=(year1-1)*10000+1215
date_3(nTimes)=(year2+1)*10000+115
i=1

do year = year1, year2
  do month = 1, 12
    if (month .eq. 2) then
      date_3(i)=year*10000+month*100+14
    else
      date_3(i)=year*10000+month*100+15
    end if
    i=i+1
  end do 
end do

; Adjust first and last date to cover full range of model run:
;date(0) = year1*10000+101
;date(i-1) = (year2+1)*10000+102

print("writing date")
a->date=date_3
print("writing datesec")
a->datesec=datesec_3
print("writing time_bnds")
a->time_bnds=time_bnds_3
print("writing CO2")
a->CO2=CO2_3
print("writing H")
a->H=H_3
print("writing NO")
a->NO=NO_3
print("writing O")
a->O=O_3
print("writing O2")
a->O2=O2_3
print("writing O2")
a->O3=O3_3
print("writing PS")
a->PS=PS_3
delete(a)
  
