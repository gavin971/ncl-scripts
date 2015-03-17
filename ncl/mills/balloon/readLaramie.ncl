load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

datapath=getenv("DATA_PATH")
rootpath=datapath+"/Deshler2014/US_Laramie_41N_105W/Nr_Stratosphere/"

print("ls "+rootpath)
filenames=systemfunc("ls "+rootpath)
;nfiles=dimsizes(filenames)
;filedates=toint(str_get_cols(filenames,0,7))
fileyears=toint(str_get_cols(filenames,0,3))

startyear=getenv("startyear")
endyear=getenv("endyear")

if (ismissing(startyear)) then
  startyear=1999
end if  

if (ismissing(endyear)) then
  endyear=2008
end if

if (startyear.gt.endyear) then
  print("FATAL: startyear = "+startyear+" > endyear = "+endyear)
  exit
end if

print("startyear = "+startyear+", endyear = "+endyear)

if (startyear.eq.endyear) then
  rangestr=startyear
else
  rangestr=startyear+"-"+endyear
end if

selectfiles=filenames(ind(fileyears.ge.startyear.and.fileyears.le.endyear))
nfiles=dimsizes(selectfiles)
outfile=datapath+"/Deshler2014/US_Laramie_41N_105W_"+rangestr+".nc"

minalt=0
maxalt=40
nalts=toint((maxalt-minalt)*2+1)
alt=fspan(minalt,maxalt,nalts)
alt!0="alt"
alt@units="km"

;bin=(/0.,0.15,0.19,0.25,0.30,0.49,0.78,0.95,1.08,1.20,1.80,2.00,2.47,4.10,4.70,9.50,15.0/)
bin=(/0.,0.15,0.25,0.30,0.49,0.78/)
bin@units="um"
bin@long_name="minimum radius"
nbins=dimsizes(bin)

;binint=new(nbins,float)
;binint(0)=0.1
;binint(1:nbins-2)=sqrt(bin(1:nbins-2)*bin(2:nbins-1))
;binint(nbins-1)=20.0

N=new((/nbins,nalts,nfiles/),"float")
N!0="bin"
N!1="alt"
N&bin=bin
N&alt=alt
N=N@_FillValue

do f = 0,nfiles-1
;do f = 0,0
  print("=======================================================================")
  f1=f+1
  print("file "+f1+" of "+nfiles)
  print(rootpath+selectfiles(f))
  header=readAsciiHead(rootpath+selectfiles(f), "  K   km")
  print((/header/))
  nhead=dimsizes(header)
  print("nhead="+nhead)
  
  lines  = asciiread(rootpath+selectfiles(f),-1,"string")
  print((/lines(nhead)/))
  ncols = str_fields_count(lines(nhead)," ")
  print("ncols="+ncols)

  fields=str_split(header(nhead-2)," ")
  print(fields)
  altcol=ind(fields.eq."Alt")+1
  print("altcol="+altcol)
  CNcol=ind(fields.eq."CN")+1
  print("CNcol="+CNcol)
  
  alt_f=tofloat(str_get_field(lines(nhead:),altcol," "))
  nalt_f=dimsizes(alt_f)
  a1=alt_f(0)
  a2=alt_f(nalt_f-1)
  print("a1="+a1+" a2="+a2)
  
  N(0,{a1:a2},f)=tofloat(str_get_field(lines(nhead:),CNcol," "))
  
  print((/lines(nhead-2)/))
  do c=1,ncols-CNcol
    lim = str_get_field(lines(nhead-2),CNcol+c," ")
    val = tofloat(str_get_field(lim,2,">"))
    print("c="+c+" lim="+lim+" val="+val)
    
    do ibin=0,nbins-1
      if (val.eq.bin(ibin)) then
        print("ibin "+ibin+", bin="+bin(ibin))
        N(ibin,{a1:a2},f)=tofloat(str_get_field(lines(nhead:),CNcol+c," "))
      end if
    end do
    delete(ibin)
    delete(lim)
    delete(val)
  end do
  delete(fields)
  delete(lines)
  delete(header)
  delete(alt_f)
end do

print("=======================================================================")

Nmedian=new((/nbins,nalts/),"float")
Nmedian!0="bin"
Nmedian!1="alt"
Nmedian&bin=bin
Nmedian&alt=alt
Nmedian=Nmedian@_FillValue

Nmin=Nmedian
Nmax=Nmedian
N25pct=Nmedian
N75pct=Nmedian

Nk=new((/nbins,nalts/),"integer")
Nk!0="bin"
Nk!1="alt"
Nk&bin=bin
Nk&alt=alt
Nk=0

print("calculating statistics...")
do k=0,nalts-1
  do b=0,nbins-1
    nsample=dimsizes(ind(.not.ismissing(N(b,k,:))))
    print("k="+k+" b="+b+" num="+nsample)
    if (nsample.gt.1) then
;      print(N(b,k,:))
      statN=stat_dispersion(N(b,k,:), False )
      Nmedian(b,k)=(/statN(8)/)
      Nmin(b,k)=(/statN(2)/)
      Nmax(b,k)=(/statN(14)/)
      N25pct(b,k)=(/statN(6)/)
      N75pct(b,k)=(/statN(10)/)
      Nk(b,k)=toint(statN(19))
    end if
  end do
end do
    

print("rm "+outfile)
rm=systemfunc("rm "+outfile)
print("creating "+outfile)
ncid=addfile(outfile,"c")
print("saving Nmedian...")
ncid->Nmedian=Nmedian
print("saving Nmin...")
ncid->Nmin=Nmin
print("saving Nmax...")
ncid->Nmax=Nmax
print("saving N25pct...")
ncid->N25pct=N25pct
print("saving N75pct...")
ncid->N75pct=N75pct
print("saving Nk...")
ncid->Nk=Nk
  
print("=======================================================================")

do k=0,nalts-1
  print(sprintf("%4.1f",alt(k))+" "+sprinti("%3.0i",Nk(0,k))+\
     "  |  "+sprintf("%3.2e",Nmedian(0,k))+" "+sprintf("%3.2e",Nmedian(1,k))+" "+sprintf("%3.2e",Nmedian(2,k))+\
     " "+sprintf("%3.2e",Nmedian(3,k))+" "+sprintf("%3.2e",Nmedian(4,k))+" "+sprintf("%3.2e",Nmedian(5,k)));+\
;     "  |  "+sprintf("%3.2e",Nmin(0,k))+" "+sprintf("%3.2e",Nmin(1,k))+" "+sprintf("%3.2e",Nmin(2,k))+\
;     " "+sprintf("%3.2e",Nmin(3,k))+" "+sprintf("%3.2e",Nmin(4,k))+" "+sprintf("%3.2e",Nmin(5,k))+\
;     "  |  "+sprintf("%3.2e",Nmax(0,k))+" "+sprintf("%3.2e",Nmax(1,k))+" "+sprintf("%3.2e",Nmax(2,k))+\
;     " "+sprintf("%3.2e",Nmax(3,k))+" "+sprintf("%3.2e",Nmax(4,k))+" "+sprintf("%3.2e",Nmax(5,k)))
end do

print("=======================================================================")

