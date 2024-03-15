%Get logical vector of proposals involving change in p and q
propChange_p = (pSeries_PROPOSAL ~= pSeries(1:end-1));
propChange_q = (qSeries_PROPOSAL ~= qSeries(1:end-1));

propChange_pq = propChange_p | propChange_q;

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

%Acceptance Rate for moves DECREASING p OR q
propDecrease_p = (pSeries_PROPOSAL < pSeries(1:end-1));
propDecrease_q = (qSeries_PROPOSAL < qSeries(1:end-1));
propDecrease_pq = propDecrease_p | propDecrease_q;
accepted_pqDecr = propDecrease_pq & accepted;
disp(['Moves decreasing p OR q:             ' num2str(sum(accepted_pqDecr)) ]);
disp(['Proposals decreasing p OR q:         ' num2str(sum(propDecrease_pq)) ]);
disp(['Acc. Rate Props decreasing p OR q:   ' num2str(sum(accepted_pqDecr)/sum(propDecrease_pq)) ]);
disp('==========================================================================');

%Acceptance Rate for moves INCREASING p
accepted_pIncr = propIncrease_p & accepted;
disp(['Moves only increased p:              ' num2str(sum(accepted_pIncr))]);
disp(['Proposals only increased p:          ' num2str(sum(propIncrease_p)) ]);
disp(['Acc. Rate only increased p:          ' num2str(sum(accepted_pIncr)/sum(propIncrease_p)) ]);
disp('==========================================================================');

%Acceptance Rate for moves INCREASING q
accepted_qIncr = propIncrease_q & accepted;
disp(['Moves only increased q:              ' num2str(sum(accepted_qIncr))]);
disp(['Proposals only increased q:          ' num2str(sum(propIncrease_q)) ]);
disp(['Acc. Rate only increased q:          ' num2str(sum(accepted_qIncr)/sum(propIncrease_q)) ]);
disp('==========================================================================');

%Acceptance Rate for moves INCREASING p AND q
propIncr_pANDq = propIncrease_p & propIncrease_q;
accepted_pqIncrAND = propIncr_pANDq & accepted;
disp(['Moves only increased p AND q:        ' num2str(sum(accepted_pqIncrAND))]);
disp(['Proposals only increased p AND q:    ' num2str(sum(propIncr_pANDq)) ]);
disp(['Acc. Rate only increased p AND q:    ' num2str(sum(accepted_pqIncrAND)/sum(propIncr_pANDq)) ]);
disp('==========================================================================');

%Acceptance Rate for moves DECREASING p
accepted_pDecr = propDecrease_p & accepted;
disp(['Moves only decreased p:              ' num2str(sum(accepted_pDecr))]);
disp(['Proposals only decreased p:          ' num2str(sum(propDecrease_p)) ]);
disp(['Acc. Rate only decreased p:          ' num2str(sum(accepted_pDecr)/sum(propDecrease_p)) ]);
disp('==========================================================================');

%Acceptance Rate for moves DECREASING q
accepted_qDecr = propDecrease_q & accepted;
disp(['Moves only decreased q:              ' num2str(sum(accepted_qDecr))]);
disp(['Proposals only decreased q:          ' num2str(sum(propDecrease_q)) ]);
disp(['Acc. Rate only decreased q:          ' num2str(sum(accepted_qDecr)/sum(propDecrease_q)) ]);
disp('==========================================================================');

%Acceptance Rate for moves DECREASING p AND q
propDecr_pANDq = propDecrease_p & propDecrease_q;
accepted_pqDecrAND = propDecr_pANDq & accepted;
disp(['Moves only decreased p AND q:        ' num2str(sum(accepted_pqDecrAND))]);
disp(['Proposals only decreased p AND q:    ' num2str(sum(propDecr_pANDq)) ]);
disp(['Acc. Rate only decreased p AND q:    ' num2str(sum(accepted_pqDecrAND)/sum(propDecr_pANDq)) ]);
disp('==========================================================================');

%Acceptance Rate for moves INCREASING p OR q while NOT DECREASING THE OTHER
propIncrease_pORq = (propIncrease_p & ~propDecrease_q) | (~propDecrease_p & propIncrease_q );
accepted_pqIncrOR = propIncrease_pORq & accepted;
disp(['Moves only increased or equal pq:    ' num2str(sum(accepted_pqIncrOR))]);
disp(['Proposals only increased or equal pq:' num2str(sum(propIncrease_pORq)) ]);
disp(['Acc. Rate only increased or equal pq:' num2str(sum(accepted_pqIncrOR)/sum(propIncrease_pORq)) ]);
disp('==========================================================================');

%Acceptance Rate for moves DECREASING p OR q while NOT INCREASING THE OTHER
propDecrease_pORq = (propDecrease_p & ~propIncrease_q) | (~propIncrease_p & propDecrease_q );
accepted_pqDecrOR = propDecrease_pORq & accepted;
disp(['Moves only decreased or equal pq:    ' num2str(sum(accepted_pqDecrOR))]);
disp(['Proposals only decreased or equal pq:' num2str(sum(propDecrease_pORq)) ]);
disp(['Acc. Rate only decreased or equal pq:' num2str(sum(accepted_pqDecrOR)/sum(propDecrease_pORq)) ]);
disp('==========================================================================');
