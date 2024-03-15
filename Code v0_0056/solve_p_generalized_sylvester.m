function [Phi_0]=solve_p_generalized_sylvester(model_info,parameters,A_tilde)

% %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% solve_p_generalized_sylvester
% inputs
% 
% 
% outputs
% 
% 
% created 10/21/2013
% 
% %%%%%%%%%%%%%%%%%%%%%%%

[Q,U] = schur(-A_tilde,'complex');
D_tilde=Q'*((parameters.X.A*parameters.X.Lambda+parameters.X.B)\parameters.X.D);
U_powers=zeros(model_info.num_endo,model_info.num_endo*(parameters.Z.p));
Phi_0=zeros(model_info.num_endo,model_info.num_exog);
for j=1:parameters.Z.p
    U_powers(:,1+(j-1)*model_info.num_endo:(j)*model_info.num_endo)=U^j;
end
x_gamma=zeros(model_info.num_endo, (model_info.num_exog)*parameters.Z.p); A=eye(model_info.num_exog,model_info.num_exog);
for j=1:parameters.Z.p
    A=A+U_powers(model_info.num_endo,model_info.num_endo+(j-1)*model_info.num_exog)*(-parameters.Z.P(:,1+(j-1)*model_info.num_exog:j*model_info.num_exog));
end
Phi_0(model_info.num_endo,:)=(-D_tilde(model_info.num_endo,:))/A;
x_gamma(model_info.num_endo,:)=Phi_0(model_info.num_endo,:)*parameters.Z.P;
for i=2:model_info.num_endo
    B=zeros(1,model_info.num_exog);A=eye(model_info.num_exog,model_info.num_exog);%Phi_0(model_info.num_endo+1-i+(j-1)*model_info.num_exog+1:model_info.num_endo+1-i+(j)*model_info.num_exog,:);
    for j=1:parameters.Z.p
        A=A+U_powers(model_info.num_endo+1-i,model_info.num_endo+1-i+(j-1)*model_info.num_endo)*(-parameters.Z.P(:,1+(j-1)*model_info.num_exog:j*model_info.num_exog));
        B=B+U_powers(model_info.num_endo+1-i,model_info.num_endo+1-i+(j-1)*model_info.num_endo+1:model_info.num_endo+(j-1)*model_info.num_endo)*x_gamma(model_info.num_endo+1-i+1:model_info.num_endo,(j-1)*model_info.num_exog+1:(j)*model_info.num_exog);
    end
%     J_index=1:parameters.Z.p;
%     A_index=model_info.num_endo+1-i+1:model_info.num_endo;
%     
% mm=A_index; 
% l_mm=length(mm); 
%       nn1=J_index(ones(l_mm,1),:);
%       nn1=nn1(:);
%        imm1=(1:l_mm);
%        imm1 = imm1(ones(parameters.Z.p,1),:);
%        mm1=mm(1,imm1.');%mm1=repmat(mm,[l_nn 1]);
%      B_index=mm1(:)+(nn1-1).*model_info.num_endo;
% %B_index=bsxfun(@plus,A_index',(J_index-1)*model_info.num_endo);B_index=B_index(:)';
% B_tilde=x_gamma(A_index,bsxfun(@plus,(J_index'-1)*model_info.num_exog,[1:model_info.num_exog]));
% B_tilde=reshape( B_tilde,[length(A_index)*parameters.Z.p, model_info.num_exog]);
% B_tilde=U_powers(model_info.num_endo+1-i,B_index)*B_tilde;
% %max(max(abs(B-B_tilde) ))  
% % max(max(abs(B)))
    Phi_0(model_info.num_endo-i+1,:)=(B-D_tilde(model_info.num_endo-i+1,:))/A;
    x_gamma(model_info.num_endo-i+1,:)=Phi_0(model_info.num_endo-i+1,:)*parameters.Z.P;
end
Phi_0=Q*Phi_0;