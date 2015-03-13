load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

datapath=getenv("DATA_PATH")
rootpath=datapath+"/Deshler2014/US_Laramie_41N_105W/Nr_Stratosphere/"

print("ls "+rootpath)
filenames=systemfunc("ls "+rootpath)
nfiles=dimsizes(filenames)

nalts=80
alt=fspan(0.5,40,nalts)
alt!0="alt"
alt@units="km"

Nave=new((/9,nalts/),"float")
Nave!0="bin"
Nave!1="alt"
Nave&alt=alt
Nave=0.0

Nsig=Nave

Nk=new((/9,nalts/),"integer")
Nk!0="bin"
Nk!1="alt"
Nk&alt=alt
Nk=0

do f = 0,nfiles-1
  print("=======================================================================")
  f1=f+1
  print("file "+f1+" of "+nfiles)
  print(rootpath+filenames(f))
  header=readAsciiHead(rootpath+filenames(f), "  K   km")
  print((/header/))
  nhead=dimsizes(header)
  print("nhead="+nhead)
  
  lines  = asciiread(rootpath+filenames(f),-1,"string")
  print((/lines(nhead)/))
  ncols = str_fields_count(lines(nhead)," ")
  print("ncols="+ncols)

  alt_f=tofloat(str_get_field(lines(nhead:),2," "))
  nalt_f=dimsizes(alt_f)
  a1=alt_f(0)
  a2=alt_f(nalt_f-1)

  cn=tofloat(str_get_field(lines(nhead:),5," "))
;  print(alt_f+" "+cn)
  
  Nk(0,{a1:a2})=Nk(0,{a1:a2})+1
  Nsig(0,{a1:a2})=Nsig(0,{a1:a2})+(Nk(0,{a1:a2})-1)/Nk(0,{a1:a2})*(cn-Nave(0,{a1:a2}))^2.0
  Nave(0,{a1:a2})=Nave(0,{a1:a2})+(cn-Nave(0,{a1:a2}))/Nk(0,{a1:a2})
  
  print((/lines(nhead-2)/))
  do c=0,ncols-6
    lim = str_get_field(lines(nhead-2),6+c," ")
    print("c="+c+" lim="+lim)
    if (lim.eq."N(r>0.15") then
      bin=1
    end if
    if (lim.eq."N(r>0.19") then
      bin=2
    end if
    if (lim.eq."N(r>0.25") then
      bin=3
    end if
    if (lim.eq."N(r>0.30") then
      bin=4
    end if
    if (lim.eq."N(r>0.49") then
      bin=5
    end if
    if (lim.eq."N(r>0.78") then
      bin=6
    end if
    if (lim.eq."N(r>1.08") then
      bin=7
    end if
    if (lim.eq."N(r>2.00") then
      bin=8
    end if
    print("bin "+bin)
    
    cn=tofloat(str_get_field(lines(nhead:),6+c," "))
    Nk(bin,{a1:a2})=Nk(bin,{a1:a2})+1
    Nsig(bin,{a1:a2})=Nsig(bin,{a1:a2})+(Nk(bin,{a1:a2})-1)/Nk(bin,{a1:a2})*(cn-Nave(bin,{a1:a2}))^2.0
    Nave(bin,{a1:a2})=Nave(bin,{a1:a2})+(cn-Nave(bin,{a1:a2}))/Nk(bin,{a1:a2})
  end do
  delete(lines)
  delete(header)
  delete(alt_f)
  delete(cn)
end do
