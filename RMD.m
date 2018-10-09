%Random matrix classification algorithm
%Alpha Lee 26/11/17
%
%Inputs
%
%training_binding: The training set for the binders 
%verification_binding: The test set for the binders 
%training_decoy: The training set for the decoys
%verification_decoy: The test set for the decoys 
%
%thres: a parameter that needs to be tuned such that the entire AUC curve
%is plotted (typically 100) 
%
%The datasets should be formatted as Nxp matrices, where N is the number of
%samples in the set and p is the number of descriptors per sample 

function AUC = RMD(training_binding,verification_binding,training_decoy,verification_decoy,thres) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%processing the binding training set 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tset_cleaned = training_binding; 

%compute z score 
[tset_cleaned_z, mu, sigma] = zscore(tset_cleaned); 
indzero_bind = find(sigma==0); %get rid of descriptors that have the same value for every member of the dataset  
tset_cleaned_z(:,indzero_bind) = []; 
mu(indzero_bind) = [];
sigma(indzero_bind) = []; 

%get covarience matrix and eigenvalues 
covar =tset_cleaned_z'*tset_cleaned_z/size(training_binding,1); 
[v, d] =eig(covar);

%Use the MP bound to get the number of significant eigenvalues  
p = size(tset_cleaned,2); 
n = size(tset_cleaned,1); 
l = diag(d);
num_eig = length(l(find(l>(1+sqrt(p/n))^2))); 

%get the significnt eigenvectors 
vv = v(:,end-num_eig+1:end); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%processing the decoy training set 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get rid of columns with the same entry 
tset_d_cleaned = training_decoy; 

%compute z score 
[tset_d_cleaned_z, mu_d, sigma_d] = zscore(tset_d_cleaned);  

indzero_bind_d = find(sigma_d==0); 
tset_d_cleaned_z(:,indzero_bind_d) = []; 
mu_d(indzero_bind_d) = []; 
sigma_d(indzero_bind_d) = []; 

%get covarience matrix and eigenvalues 
covar_d =tset_d_cleaned_z'*tset_d_cleaned_z/size(training_decoy,1); 
[v_decoy, d_decoy] =eig(covar_d);

%Use the MP bound to get the number of significant eigenvalues  
p = size(tset_d_cleaned,2); 
n = size(tset_d_cleaned,1); 
l_decoy = diag(d_decoy);
num_eig_d = length(l_decoy(find(l_decoy>(1+sqrt(p/n))^2))); 

%get the significnt eigenvectors 
vv_decoy = v_decoy(:,end-num_eig_d+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%processing the binding verification set 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first look at how close is the verification set to the binding set
verification_binding1 = verification_binding; 
verification_binding2 = verification_binding; 
verification_binding1(:,indzero_bind) = [];
verification_binding2(:,indzero_bind_d) = [];

%first look at how close the compounds are to the active training set 

%mean center and scale the verification set w.r.t. the active training set
veriset_mu = (verification_binding1-repmat(mu,size(verification_binding1,1),1))./repmat(sigma,size(verification_binding1,1),1);   
coeff = veriset_mu*vv; 

%project back into the ligand space 
proj_vect = (vv*coeff')';
norm_test = sqrt(sum((proj_vect-veriset_mu).^2,2));

%now look at how close the compounds are to the decoy training set 

%mean center and scale the verification set w.r.t. the decoy training
%set
veriset_mu = (verification_binding2-repmat(mu_d,size(verification_binding2,1),1))./repmat(sigma_d,size(verification_binding2,1),1);   
coeff = veriset_mu*vv_decoy; 

%project back into the ligand space 
proj_vect = (vv_decoy*coeff')';
norm_test_neg = sqrt(sum((proj_vect-veriset_mu).^2,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%processing the decoy verification set 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verification_decoy1 = verification_decoy;
verification_decoy2 = verification_decoy;
verification_decoy1(:,indzero_bind) = [];
verification_decoy2(:,indzero_bind_d) = [];

%first look at how close the compounds are to the active training set 

%mean center and scale the verification set w.r.t. the active training set
veriset_d_mu = (verification_decoy1-repmat(mu,size(verification_decoy1,1),1))./repmat(sigma,size(verification_decoy1,1),1);   
coeff_d = veriset_d_mu*vv; 

%project back into the ligand space 
proj_vect_decoy = (vv*coeff_d')';
norm_test_decoy = sqrt(sum((proj_vect_decoy-veriset_d_mu).^2,2));

%now look at how close the compounds are to the decoy training set 

%mean center and scale the verification set w.r.t. the decoy training
%set
veriset_d_mu = (verification_decoy2-repmat(mu_d,size(verification_decoy2,1),1))./repmat(sigma_d,size(verification_decoy2,1),1);   
coeff_d = veriset_d_mu*vv_decoy; 

%project back into the ligand space 
proj_vect_decoy = (vv_decoy*coeff_d')';
norm_test_decoy_neg = sqrt(sum((proj_vect_decoy-veriset_d_mu).^2,2));

threshold = -thres:0.01:thres;  
for ii = 1:length(threshold) 
% compute false negative and false positve 

     true_pos(ii) = length(find(norm_test < (norm_test_neg + threshold(ii))))/length(norm_test);  
     false_pos(ii) = length(find(norm_test_decoy < (norm_test_decoy_neg + threshold(ii)) ))/ length(norm_test_decoy); 

end

AUC = trapz(false_pos,true_pos);

end  