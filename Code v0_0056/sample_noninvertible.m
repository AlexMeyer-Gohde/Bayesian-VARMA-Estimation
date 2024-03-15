function [parameters]=sample_noninvertible(parameters)

%%%%%%%%%%%%%%%%%%%%%%
%  Only works for now with scalar series
%  There can be multiple series, but the Matrices in parameters.Z.Q must be
%  diagonal
%  Additionallay, the innovations can be correlated though
%%%%%%%%%%%%%%%%%%

[num_exog temp]=size(parameters.Z.Q);
q=temp/num_exog;
 Sigma_factor=zeros(1,num_exog);
% for j=0:q
%   is_diag(j+1)=is_diag(parameters.Z.Q(:,j*num_exog+1:(j+1)*num_exog));
% end
% if is_diag(is_diag==0)
% 	disp('Sorry, sample_noninvertible only works on scalar armas');
% else
for i=1:num_exog
    try
	Q=[1 parameters.Z.Q(i,i+(0:num_exog:(q-1)*num_exog))];
	[q_i]=find(Q~=0); q_i=max(q_i)-1;
    Q=Q(1:q_i+1);
    Q_roots=cplxpair(roots(Q));
    number_of_complex_roots=sum(logical(imag(Q_roots)));
    pairs_of_complex_roots=number_of_complex_roots/2;
    set=0:q_i-pairs_of_complex_roots;
    number_of_flips=randsample(set, 1, true, factorial(max(set))./(factorial(max(set)-set).*factorial(set)));
    if number_of_flips>0
        stacked_roots=[1:2:number_of_complex_roots-1, number_of_complex_roots+1:q;2:2:number_of_complex_roots,zeros(1,q-number_of_complex_roots)];
        stacked_roots=stacked_roots(:,randperm(pairs_of_complex_roots+q-number_of_complex_roots));
        which_roots_to_flip=stacked_roots(:,1:number_of_flips);
        which_roots_to_flip=which_roots_to_flip(which_roots_to_flip~=0);
        Sigma_factor(i)=abs(prod(Q_roots(which_roots_to_flip)));
        Q_roots(which_roots_to_flip)=1./Q_roots(which_roots_to_flip);
        poly_q_roots=poly(Q_roots);
        parameters.Z.Q(i,i+(0:num_exog:(q-1)*num_exog))=poly_q_roots(2:end);
    else
        which_roots_to_flip=[];
        Sigma_factor(i)=1;
    end
    catch error
        disp('Boom!')
    end
end
parameters.Z.Sigma=Sigma_factor'*Sigma_factor.*parameters.Z.Sigma;
end
