; Constants
R=8.31441 ; J/mol/k = Pa m3/mol/K
pi=acos(0.0)*2.0
sqrt2=2.^0.5
Na = 6.022e23 ; Avogadro's number (/mol)
d = 4.18e-10 ; molecular diameter in meters
MH2SO4=98.078479e-3 ; kg/mol

a=addfile("/glade/p/work/mmills/inputdata/phot/temp_prs_GT200nm_JPL10_expH2SO4_QY_c140624.nc","w")
print("reading jh2so4")
jh2so4=a->jh2so4

; Pressure levels used in the cross section table:
press = a->pressure
npress=dimsizes(press)
press = press*100. ; hPa -> Pa

; Temperatures used in the cross section table:
temp = a->temperature
ntemp=dimsizes(temp)

; statistical k_rrkm rates from Miller et al. (GRL, 2007) for v = 6 (519 nm), 
; v = 5 (606 nm), v = 4 (741 nm):
rrkm = (/ 1.1e-9, 8.9e-9, 1.7e-7 /) ; sec
nwave = dimsizes(rrkm)
rrkm = 1./rrkm ; sec-1

xs = (/1.43e-26, 1.8564e-25, 3.086999e-24 /)

do p=0,npress-1
  print("press="+press(p))
  do t=0,ntemp-1
    if (p.eq.npress-1)
      qy = (/ 1.0, 1.0, 1.0 /); full cross section at top level
    else
      lambda=R*temp(t)/(sqrt2*pi*d*d*Na*press(p))
      veloc=(8.*R*temp(t)/(pi*MH2SO4))^0.5
      fcollis=veloc/lambda ; collision rate (s-1)
      qy=rrkm(:)/(rrkm(:)+fcollis)
;      print("temp="+temp(t))
;      print("lambda="+lambda)
;      print("veloc="+veloc)
;      print("fcollis="+fcollis)
;      print("qy="+qy)
    end if
    jh2so4(p,t,62)=xs(0)*qy(0)
    jh2so4(p,t,64)=xs(1)*qy(1)
    jh2so4(p,t,66)=xs(2)*qy(2)
  end do
end do

print("writing jh2so4")
a->jh2so4 = (/ jh2so4 /)
