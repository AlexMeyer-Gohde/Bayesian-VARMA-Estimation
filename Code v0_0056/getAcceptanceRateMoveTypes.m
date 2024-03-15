%Get logical vector of proposals involving change in p and q
propChange_p = (pSeries_PROPOSAL ~= pSeries(1:end-1));
propChange_q = (qSeries_PROPOSAL ~= qSeries(1:end-1));

%Get logical vector of acceptance
accepted = (logPosteriorSeries(1:end-1) ~= logPosteriorSeries(2:end));

%Total acceptance ratio
disp(['Total Acceptance Rate:               ' num2str(sum(accepted)/length(accepted)) ]);
disp('==========================================================================');

%Acceptance Rate for within-model moves
accepted_nochange = ~propChange_pq & accepted;
disp(['Within-Model Moves:                  ' num2str(sum(accepted_nochange)) ]);
disp(['Within-Model Proposals:              ' num2str(sum(~propChange_pq)) ]);
disp(['Within-Model Move Acceptance Rate:   ' num2str(sum(accepted_nochange) / sum(~propChange_pq)) ]);
disp('==========================================================================');

%Acceptance Rate for moves changing either p or q
propChange_pq = (propChange_p | propChange_q);
accepted_pq = propChange_pq & accepted;
disp(['Model Moves:                         ' num2str(sum(accepted_pq)) ]);
disp(['Model Proposals:                     ' num2str(sum(propChange_pq)) ]);
disp(['Model Move Acceptance Rate:          ' num2str(sum(accepted_pq)/sum(propChange_pq)) ]);
disp('==========================================================================');

%Acceptance Rate for moves INCREASING p OR q (logical)
propIncrease_p = (pSeries_PROPOSAL > pSeries(1:end-1));
propIncrease_q = (qSeries_PROPOSAL > qSeries(1:end-1));
propIncrease_pq = propIncrease_p | propIncrease_q;
accepted_pqIncr = propIncrease_pq & accepted;
disp(['Moves increasing p OR q:             ' num2str(sum(accepted_pqIncr)) ]);
disp(['Proposals increasing p OR q:         ' num2str(sum(propIncrease_pq)) ]);
disp(['Acc. Rate Props increasing p OR q:   ' num2str(sum(accepted_pqIncr)/sum(propIncrease_pq)) ]);
disp('==========================================================================');

