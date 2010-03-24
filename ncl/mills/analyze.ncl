; Meteoritic Dust Analysis
;
; Add some analysis fields to history files created by the WACCM/CARMA
; meteoritic dust model.


WACCM_RUN = getenv("WACCM_RUN")

if (ismissing(WACCM_RUN)) then
  WACCM_RUN = "results"
end if
print("Using WACCM_RUN: " + WACCM_RUN)
print("")
print("Analyzing History Files")

;HISTORYPATTERN = "run/" + WACCM_RUN + "/*cam2.h0*"
HISTORYPATTERN = "run/" + WACCM_RUN + "/h0.2034.ave.nc"

load "ncl-scripts/dust.ncl"

analyzeHistory(HISTORYPATTERN)
