;load "scripts/SootTempChange.ncl"

;saveMeans("Soot5Tg", "PassiveSoot")
;plotTempO3Change("Soot5Tg", "PassiveSoot")
;plotLOx("Soot5Tg", "PassiveSoot")

load "scripts/PlotO3Depl.ncl"
plotO3Depletion("Soot5Tg", "PassiveSoot")
