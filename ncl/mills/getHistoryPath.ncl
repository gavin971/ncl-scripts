function getHistoryPath ( historyRun:string )
local historyPath
begin

  if (isfilepresent("/Volumes/ptmp")) then
    historyPath="/Volumes/ptmp/mills/waccm/dust-sulfur/run/"+historyRun+"/"
  else
    if (isfilepresent("csm")) then
      historyPath="csm/"+historyRun+"/atm/hist/"
    else
      if (isfilepresent("waccm")) then
        historyPath="waccm/dust-sulfur/run/"+historyRun+"/"
      else
        if (isfilepresent("dust-sulfur")) then
          historyPath="dust-sulfur/run/"+historyRun+"/"
        else
          if (isfilepresent("run")) then
            historyPath="run/"+historyRun+"/"
          else
            if (isfilepresent(historyRun)) then
              historyPath=historyRun+"/"
            end if
          end if
        end if
      end if
    end if
  end if

  return(historyPath)
end
