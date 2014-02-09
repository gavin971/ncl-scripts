load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

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


sigmag=(/1.6, 1.6, 1.8/) ; geometric standard deviations for 3 MAM3 modes

historyPattern=rootpath+caseid+"/atm/hist/*.cam.h0.199[23]*"
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


