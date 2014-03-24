load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("calcSizeDistN2O")
procedure calcSizeDistN2O(caseid:string, dgVar:string, yearPattern:string,\
                          yearString:string)

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

  N2Oppb = (/ 300, 250, 200, 150, 100, 50 /)
  N2Oppb@long_name="N2O boundaries"
  N2Oppb@units="ppbv"
  
  nBounds = dimsizes(N2Oppb)
  nRanges = nBounds-1
  
  dlog10Dp = 0.056268 ; bin width of Wilson FCAS measurements

  ; Bin geometric mean diameters for Wilson FCAS measurements:
  DpMean = (/ 0.064015516,0.072870752,0.082950928,0.094425491,0.107487324, \
              0.122355994,0.139281441,0.158548176,0.18048007,0.205445792, \
              0.233865011,0.266215446,0.303040902,0.344960406,0.392678615, \
              0.446997661,0.508830634,0.579216934,0.659339737,0.750545889, \
              0.854368544,0.972552939,1.107085724,1.260228365,1.434555155, \
              1.632996487,1.858888114,2.116027222,2.40873626,2.741935601, \
              3.121226247,3.5, 4. , 4.5, 5. /)
  DpMean@units="um"
  
  nBins = dimsizes(DpMean)
  log10DpMean = log10(DpMean)

;  sigmag=(/1.6, 1.6, 1.8/) ; geometric standard deviations for 3 MAM3 modes
  sigmag=(/1.6, 1.6, 1.2/) ; geometric standard deviations for 3 MAM3 modes
  log10sigmag=log10(sigmag)
    
  N2Oppb2=N2Oppb(0:nRanges-1)
  
  dVdlog10Dp_avg=new((/nRanges,nBins/),"float")
  dVdlog10Dp_avg!0="N2O"
  dVdlog10Dp_avg&N2O=N2Oppb2
  dVdlog10Dp_avg!1="diameter"
  dVdlog10Dp_avg&diameter=DpMean
  
  dVdlog10Dp_std=dVdlog10Dp_avg
  
  outfilename=rootpath+caseid+"/atm/proc/SizeDistN2O_"+yearString+".nc"
  dum=systemfunc("rm "+outfilename)
  outfile=addfile(outfilename,"c")

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
  N2Obounds@units=N2O@units
  N2Obounds@long_name="N2O boundaries"
  
  do n = 0, nRanges-1
    np1=n+1
    print("    ")
    print("=========================================================================")
;    print(N2Oppb)
;    print(n+" "+np1)
    print(np1+" of "+nRanges+": N2O range "+N2Oppb(n)+" to "+N2Oppb(np1))
    if (N2Oppb(n).le.N2Oppb(n+1)) then
      print("calcSizeDistN2O error: N2Oppb must be in descending order. Aborting.")
      print(N2Oppb)
      exit
    end if
    
    iN2O = ind(N2O_1D.lt.N2Obounds(n).and.N2O_1D.gt.N2Obounds(n+1))

;    printVarSummary(iN2O)
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
  
  print("writing output")
  outfile->N2Oppb=N2Oppb
  outfile->N2Obounds=N2Obounds
  outfile->dVdlog10Dp_avg=dVdlog10Dp_avg
  outfile->dVdlog10Dp_std=dVdlog10Dp_std
  
end

;caseid = "b.e12.BRCP45W5CN.f19_g16.MAMstratchem.003"
;yearPattern = "200[67]"
;yearString = "2006-07"

;caseid = "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-25km.15S-15N.zavg.QBOon.003"
;caseid = "b.e12.B55TRW5CN.f19_g16.Pin15Tg.16-20km.QBOon.001"

;caseid = "b.e12.B55TRW5CN.f19_g16.Pin20Tg.22-26km.1col.QBOon.001"
;caseid = "f.e12.FSDW5.f19_f19.Pin15Tg.20-27km.0-18N.96-118E.QBOon.001"
;caseid = "f.e12.FSDW5.f19_f19.Pin10Tg.20-25km.15S-15N.zavg.QBOon.002"
;caseid = "f.e12.FSDW5.f19_f19.Pin20Tg.19-27km.0-20N.120E.QBOon.001"
;caseid = "f.e12.FSDW5.f19_f19.Pin10Tg.20-25km.15S-15N.zavg.a3dgnl44.001"
caseid = "f.e12.FSDW5.f19_f19.Pin15Tg.19-27km.0N-20N.120E.a3dgnl44.001"
;yearPattern = "199[23]"
;yearString = "1992-93"
yearPattern = "199[67]"
yearString = "1996-97"

;dgVar = "dgnumwet"
dgVar = "dgnd_a0"
calcSizeDistN2O(caseid, dgVar, yearPattern, yearString)

