load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("calcSizeDistN2O")
procedure calcSizeDistN2O(caseid:string, dgVar:string, yearPattern:string,\
                          yearString:string,N2Oppb[*]:integer, \
                          DpMean[*]:float, dVdlog10Dp_avg[*,*]:float, \
                          dVdlog10Dp_std[*,*]:float)

; input variables
;    caseid      : e.g "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-25km.15S-15N.zavg.QBOon.003"  
;    dgVar       : "dgnumwet" or "dgnd_a0"
;    yearPattern : e.g. "199[23]"
;    yearString  : e.g. "1992-93"
;    N2Oppb      : N2O bounds in descending order, e.g. (/ 300, 250, 200, 150, 100, 50 /)

; output variables
;    DpMean         : Bin geometric mean diameters
;    dVdlog10Dp_avg : normalized mean volume distribution 
;    dVdlog10Dp_std : normalized standard deviation about the mean volume distribution
  
begin
  pi      = 4.0 * atan(1.0) 
  sqrt2pi=sqrt(2*pi)

  rootpath = "/glade/scratch/mmills/archive/"

  nBounds = dimsizes(N2Oppb)
  
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
  
  dVdlog10Dp_avg=new((/nBounds-1,nBins/),"float")
  dVdlog10Dp_std=new((/nBounds-1,nBins/),"float")

  historyPattern=rootpath+caseid+"/atm/proc/h0zm/*.cam.h0zm."+yearPattern+"*"
  historyNames = systemfunc("ls " + historyPattern)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  print("historyNames:"+dimsizes(historyNames))
  print("reading N2O")
  N2O = addfiles_GetVar(history, historyNames, "N2O")

  ; dry diameter not saved, but wet should be the same as dry for these runs
  print("reading "+dgVar+"1")
  dgnum1 = addfiles_GetVar(history, historyNames, dgVar+"1")
  print("reading "+dgVar+"2")
  dgnum2 = addfiles_GetVar(history, historyNames, dgVar+"2")
  print("reading "+dgVar+"3")
  dgnum3 = addfiles_GetVar(history, historyNames, dgVar+"3")
  print("reading num_a1")
  num_a1 = addfiles_GetVar(history, historyNames, "num_a1")
  print("reading num_a2")
  num_a2 = addfiles_GetVar(history, historyNames, "num_a2")
  print("reading num_a3")
  num_a3 = addfiles_GetVar(history, historyNames, "num_a3")

  dg1=ndtooned(dgnum1)
  dg2=ndtooned(dgnum2)
  dg3=ndtooned(dgnum3)

  n1=ndtooned(num_a1)
  n2=ndtooned(num_a2)
  n3=ndtooned(num_a3)

  printVarSummary(N2O)
  ds=dimsizes(N2O)

  N2O_1D = ndtooned(N2O)
  N2Obounds=N2Oppb*1.e-9
  
  do n = 0, nBounds-1
  
    print(n+" of "+nBounds+": N20 range "+N2Oppb[n]+" to "+N2Oppb[n+1])
    if (N2Oppb[n].le.N2Oppb[n+1]) then
      print("calcSizeDistN2O error: N2Oppb must be in descending order. Aborting.")
      print(N2Oppb)
      exit
    end if
    
    iN2O = ind(N2O_1D.lt.N2Obounds[n].and.N2O_1D.gt.N2Obounds[n+1])

    printVarSummary(iN2O)
    nN2O = dimsizes(iN2O)

    log10dgvol=new((/3,nN2O/),"float")
    numarr=log10dgvol

    log10dgvol(0,:) = log10(dg1(iN2O)*1e6) + 3*log10sigmag(0)*log10sigmag(0)
    log10dgvol(1,:) = log10(dg2(iN2O)*1e6) + 3*log10sigmag(1)*log10sigmag(1)
    log10dgvol(2,:) = log10(dg3(iN2O)*1e6) + 3*log10sigmag(2)*log10sigmag(2)

    numarr(0,:) = ndtooned(n1(iN2O))
    numarr(1,:) = ndtooned(n2(iN2O))
    numarr(2,:) = ndtooned(n3(iN2O))

    sizedist=new((/nBins,nN2O/),"float")
    sizedist(:,:)=0.0

    do t = 0, nN2O-1
      if (mod(t,1000).eq.0) then
        print(t+" of "+nN2O)
      end if
      do i = 0, nBins-1
        do m=0,2 ; modes
          dVdlog10Dp=3.0*log10DpMean(i)+(9./2.)*log10sigmag(m)*log10sigmag(m)
          dVdlog10Dp=dVdlog10Dp-(((log10DpMean(i)-log10dgvol(m,t))/log10sigmag(m))^2.)/2
          dVdlog10Dp=numarr(m,t)/DpMean(i)/log10sigmag(m)*10^dVdlog10Dp
          sizedist(i,t)=sizedist(i,t)+dVdlog10Dp
        end do
      end do
      ; normalize distribution
      sizedist(:,t)=sizedist(:,t)/sum(sizedist(:,t))/dlog10Dp
    end do

    print("calculating average")
    dVdlog10Dp_avg(n,:)=dim_avg(sizedist)
    print("calculating standard deviation")
    dVdlog10Dp_std(n,:)=dim_stddev(sizedist)
    
    delete(iN2O)
    delete(nN2O)
    delete(log10dgvol)
    delete(numarr)
    delete(sizedist)
    delete(dVdlog10Dp)
    
  end do
end

undef("saveSizeDistN2O")
procedure saveSizeDistN2O
  
end
