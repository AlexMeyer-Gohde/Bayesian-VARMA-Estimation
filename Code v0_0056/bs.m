cntrCorrect = 0;
for cntr = 1 : length(consolidatedResults)
   if consolidatedResults(cntr).pPostMax == 3 && consolidatedResults(cntr).qPostMax == 2
       cntrCorrect = cntrCorrect +1;
   end;
   disp(cntrCorrect/length(consolidatedResults));
end