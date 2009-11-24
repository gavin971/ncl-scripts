load "tests/common2.ncl"
load "tests/getHistoryPath.ncl"

targetLat = 78
targetLon = 0
targetLev=0.00230591145   ; 86.7 km at 78N
targetDate = "19950615"

print("Lat="+targetLat)
print("Lev="+targetLev)

;historyRun="125tabaz22"
;historyPath = getHistoryPath(historyRun) + "alldaily/"

;file1="all."+historyRun+".cam2.h3.1995-06-14-00000.nc"
;history=addfile(historyPath+file1,"r")

 ; vars(time, lev, lat, lon, bins)
;  print("  reading MIXSULNDBIN")
;  mixsulndbin  = history->MIXSULNDBIN
;  print("  reading SULFRAD")
;  sulfrad  = history->SULFRAD

;  print("  reading PURSULNDBIN")
;  pursulndbin  = history->PURSULNDBIN
;  print("  reading PURERAD")
;  purerad  = history->PURERAD
;  print("  reading SULDEN")
;  sulden  = history->SULDEN
;  t  = history->T
;  wtpct=history->WTPCT
;  lev = history->lev
;  lat = history->lat
;  lon = history->lon
;  time  = history->time
;  date  = history->date

;  iTime=ind(date.eq.targetDate)

;  rad1 = purerad(iTime,{targetLev},{targetLat},{targetLon},:)
;  rad2 = rad1
;  rhoaer = sulden(iTime,{targetLev},{targetLat},{targetLon})
;  temp = t(iTime,{targetLev},{targetLat},{targetLon})
;  wp = wtpct(iTime,{targetLev},{targetLat},{targetLon})
 
temp= 149.330815851386916
wp=45.6412523721937973
rhoaer=1.572573157
rad1=(/1.08E-08,1.36E-08,1.72E-08,2.17E-08,2.73E-08,3.44E-08/)
 
  
  dnwtp = (/0,1,5,10,20,25,30,35,40,41,45,50,53,55,56,60,65,66,70, \
   72,73,74,75,76,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93, \
   94,95,96,97,98,100/)
  dnc0 = (/1,1.13185,1.17171,1.22164,1.3219,1.37209,1.42185,1.4705, \
   1.51767,1.52731,1.56584,1.61834,1.65191,1.6752,1.68708,1.7356,   \
   1.7997,1.81271,1.86696,1.89491,1.9092,1.92395,1.93904,1.95438,   \
   1.98574,2.00151,2.01703,2.03234,2.04716,2.06082,2.07363,2.08461, \
   2.09386,2.10143,2.10764,2.11283,2.11671,2.11938,2.12125,2.1219,  \
   2.12723,2.12654,2.12621,2.12561,2.12494,2.12093/)
  dnc1 = (/0,-0.000435022,-0.000479481,-0.000531558,-0.000622448,  \
   -0.000660866,-0.000693492,-0.000718251,-0.000732869,-0.000735755,\
   -0.000744294,-0.000761493,-0.000774238,-0.00078392,-0.000788939, \
   -0.00080946,-0.000839848,-0.000845825,-0.000874337,-0.000890074, \
   -0.00089873,-0.000908778,-0.000920012,-0.000932184,-0.000959514, \
   -0.000974043,-0.000988264,-0.00100258,-0.00101634,-0.00102762,   \
   -0.00103757,-0.00104337,-0.00104563,-0.00104458,-0.00104144,     \
   -0.00103719,-0.00103089,-0.00102262,-0.00101355,-0.00100249,     \
   -0.00100934,-0.000998299,-0.000990961,-0.000985845,-0.000984529, \
   -0.000989315/)
   
  i = 0  
  print("")
  print("weight percent="+wp)
  do while (wp .gt. dnwtp(i))
    i = i + 1
    print("i="+i+" dnwtp(i)="+dnwtp(i))
  end do
  print("i="+i)
  den1=dnc0(i-1)+dnc1(i-1)*temp
  den2=dnc0(i)+dnc1(i)*temp
  drho_dwt=(den2-den1)/(dnwtp(i)-dnwtp(i-1))

  if (i.eq.1 .or. wp.eq.dnwtp(i)) then
    rho=den2
  else
    frac=(dnwtp(i)-wp)/(dnwtp(i)-dnwtp(i-1))
    rho=den1*frac+den2*(1.0-frac)
  end if
  
  print("dnwtp(i)="+dnwtp(i))
  print("dnwtp(i+1)="+dnwtp(i+1))
  print("rho [g/cm2]="+rho)
  print("rhoaer [g/cm3]="+rhoaer)  
  print("drho_dwt ="+drho_dwt)
  print("")
  
  stwtp = (/0, 23.8141, 38.0279, 40.6856, 45.335, 52.9305, 56.2735, \
      59.8557, 66.2364, 73.103, 79.432, 85.9195, 91.7444, 97.6687, 100/)

  stc0 = (/117.564, 103.303, 101.796, 100.42, 98.4993, 91.8866,    \
      88.3033, 86.5546, 84.471, 81.2939, 79.3556, 75.608, 70.0777, \
      63.7412, 61.4591 /)

  stc1 = (/-0.153641, -0.0982007, -0.0872379, -0.0818509, \
      -0.0746702, -0.0522399, -0.0407773, -0.0357946, -0.0317062,  \
      -0.025825, -0.0267212, -0.0269204, -0.0276187, -0.0302094,   \
      -0.0303081 /)
      
  i = 0  
  print("")
  print("weight percent="+wp)
  do while (wp .gt. stwtp(i))
    i = i + 1
    print("i="+i+" stwtp(i)="+stwtp(i))
  end do
  print("i="+i)
  
  sig1=stc0(i-1)+stc1(i-1)*temp
  sig2=stc0(i)+stc1(i)*temp
  dsig_dwt=(sig2-sig1)/(stwtp(i)-stwtp(i-1))

  if (i.eq.0 .or. wp.eq.stwtp(i)) then
    surf_tens=sig2
  else
    frac=(stwtp(i)-wp)/(stwtp(i)-stwtp(i-1))
    surf_tens=sig1*frac+sig2*(1.0-frac)
  end if
  
  print("stwtp(i)="+stwtp(i))
  print("stwtp(i-1)="+stwtp(i-1))
  print("surf tens [erg/cm2]="+surf_tens)
  print("dsig_dwt="+dsig_dwt)
  print("")
  print("temp="+temp)
  print("MWH2O [g/mole]="+MWH2O)
  print("RGAS [J/K/kmole]="+RGAS)
  RGAS2 = RGAS * 1e7 / 1e3
  print("RGAS2 [erg/K/mole]="+RGAS2)

  kelvinH2O_a = 2 * MWH2O * surf_tens / (rho * RGAS2 * temp * rad1)
  kelvinH2O_b = 1+ wp*drho_dwt/rho - 3.*wp*dsig_dwt/(2.*surf_tens)
  print("")
  print("kelvinH2O_a="+kelvinH2O_a)
  print("kelvinH2O_b="+kelvinH2O_b)
  print("")
  
  kelvinH2O_1 = exp (kelvinH2O_a)
  kelvinH2O = exp (kelvinH2O_a*kelvinH2O_b)
  
  print("radius (nm)  |  preexponential  |  Kelvin factor 1  |  Kelvin factor")
  do i = 0, 37
    print(rad1(i)*1e7+" | "+kelvinH2O_a(i)+" | "+kelvinH2O_1(i)+" | "+kelvinH2O(i))
  end do
