load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

pi      = 4.0 * atan(1.0) 
sqrt2pi=sqrt(2*pi)

caseid = "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-25km.15S-15N.zavg.QBOon.002"
rootpath = "/glade/scratch/mmills/archive/"

dlog10Dp = 0.056268 ; bin width of Wilson FCAS measurements

; Bin geometric mean diameters for Wilson FCAS measurements:
DpMean = (/ 0.064015516,0.072870752,0.082950928,0.094425491,0.107487324, \
            0.122355994,0.139281441,0.158548176,0.18048007,0.205445792, \
            0.233865011,0.266215446,0.303040902,0.344960406,0.392678615, \
            0.446997661,0.508830634,0.579216934,0.659339737,0.750545889, \
            0.854368544,0.972552939,1.107085724,1.260228365,1.434555155, \
            1.632996487,1.858888114,2.116027222,2.40873626,2.741935601, \
            3.121226247 /)
nBins = dimsizes(DpMean)
log10DpMean = log10(DpMean)

sigmag=(/1.6, 1.6, 1.8/) ; geometric standard deviations for 3 MAM3 modes
log10sigmag=log10(sigmag)

historyPattern=rootpath+caseid+"/atm/proc/h0zm/*.cam.h0zm.199[23]*"
historyNames = systemfunc("ls " + historyPattern)
history = addfiles(historyNames, "r")
ListSetType(history, "cat")
print("historyNames:"+dimsizes(historyNames))
print("reading N2O")
N2O = addfiles_GetVar(history, historyNames, "N2O")

; dry diameter not saved, but wet should be the same as dry for these runs
print("reading dgnumwet1")
dgnum1 = addfiles_GetVar(history, historyNames, "dgnumwet1")
print("reading dgnumwet2")
dgnum2 = addfiles_GetVar(history, historyNames, "dgnumwet2")
print("reading dgnumwet3")
dgnum3 = addfiles_GetVar(history, historyNames, "dgnumwet3")
print("reading num_a1")
num_a1 = addfiles_GetVar(history, historyNames, "num_a1")
print("reading num_a2")
num_a2 = addfiles_GetVar(history, historyNames, "num_a2")
print("reading num_a3")
num_a3 = addfiles_GetVar(history, historyNames, "num_a3")

printVarSummary(N2O)
ds=dimsizes(N2O)

N2O_1D = ndtooned(N2O)
i250 = ind(N2O_1D.lt.300.e-9.and.N2O_1D.gt.250.e-9)

printVarSummary(i250)
n250 = dimsizes(i250)

log10dgvol=new((/3,n250/),"float")
num=dgnum

log10dgvol(0,:) = ndtooned(log10(dgnum1*1e6)) + 3*log10sigmag(0)*log10sigmag(0)
log10dgvol(1,:) = ndtooned(log10(dgnum2*1e6)) + 3*log10sigmag(1)*log10sigmag(1)
log10dgvol(2,:) = ndtooned(log10(dgnum3*1e6)) + 3*log10sigmag(2)*log10sigmag(2)
num(0,:) = ndtooned(num_a1)
num(1,:) = ndtooned(num_a2)
num(2,:) = ndtooned(num_a3)

sizedist=new((/nBins,n250/),"float")
sizedist(:,:)=0.0

do t = 0, n250-1
  print(t+" of "+n250)
  do i = 0, nBins-1
    do m=0,2 ; modes
      dVdlog10Dp=3.0*log10DpMean(i)+(9./2.)*log10sigmag(m)*log10sigmag(m)
      dVdlog10Dp=dVdlog10Dp-(((log10DpMean(i)-log10dgvol(m,t))/log10sigmag(m))^2.)/2
      dVdlog10Dp=num(m,t)/DpMean(i)/log10sigmag(m)*10^dVdlog10Dp
      sizedist(i,t)=sizedist(i,t)+dVdlog10Dp
    end do
  end do
end do
sizedist(i,t)=sizedist(i,t)/sqrt2pi

print("calculating average")
dVdlog10Dp_avg=dim_avg(sizedist)
print("calculating standard deviation")
dVdlog10Dp_std=dim_stddev(sizedist)

pdffile="SizeDistN2O"
wks   = gsn_open_wks ("ps", "$IMAGE_PATH/"+pdffile)                   ; open workstation
res=True
res@xlog=True

plot = gsn_csm_xy (wks, DpMean, dVdlog10Dp_avg, res)

