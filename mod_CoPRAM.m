%% CoPRAM
%"Fast, sample-efficient algorithms for structured phase retrieval"
% Gauri Jagatap and Chinmay Hegde, Iowa State University
% email: gauri@iastate.edu

% recover x given y_abs, A.
% y_abs = |Ax|.
% A is sensing matrix (m x n); Normal(0,1) distributed; typically m < n.
% x is s-sparse (n x 1).
% y_abs are magnitude only measurements (m x 1).
function [x,err_hist,p] =  mod_CoPRAM(y_mod,A,s,max_iter,R,tol1,tol2,z,noise)
%%updated 5/31/2017

% %% initialize parameters
 [m,n] = size(A);
% %If ground truth is unknown
% if nargin < 9
%     z = zeros(n,1);
% end
% p = zeros(m,1); %phase vector
% error_hist(1,1) = 1; %stores error in measurement model
% error_hist(1,2) = 1; %stores relative error b/w iterations
% Marg = zeros(1,n); %marginals
% MShat = zeros(s); %truncated correlation matrix
% AShat = zeros(m,s); %truncated sensing matrix
% supp = zeros(1,n); %indicator for initial support Shat
% y_mod2 = y_mod.^2; %quadratic measurements
% phi_sq = sum(y_mod2)/m;
% phi = sqrt(phi_sq); %signal power
% 
% %% s-Truncated sensing vectors
% 
% %signal marginals
% Marg = ((y_mod2)'*(A.^2))'/m; % n x 1
% [Mg MgS] = sort(Marg,'descend');
% S0 = MgS(1:s); %pick top s-marginals
% Shat = sort(S0); %store indices in sorted order
% %supp(Shat) = 1; figure; plot(supp); %support indicator
% AShat = A(:,Shat); % m x s %sensing sub-matrix
% 
% %% Truncated measurements
% %TAF = 'on'; %consider only truncated measurements m' < m ; gives marginally better performance 
% TAF = 'off'; %consider all measurements = m ; aligns with code presented in paper
% 
% switch TAF
%     case 'on'
%         card_Marg = ceil(m/6);
%         %large measurements - amplitude flow
%         for i=1:m
%             M_eval(i) = y_mod(i)/norm(AShat(i,:));
%         end 
%         [~,MmS] = sort(M_eval,'descend');
%         Io = MmS(1:card_Marg); %indices between 1 to m
%     case 'off'
%         card_Marg = m;
%         Io = 1:card_Marg;
% end
% 
% %% initialize x
% %compute top singular vector according to thresholded sensing vectors and large measurements
% for i = 1:card_Marg
%     ii = Io(i);
%     MShat = MShat + (y_mod2(ii))*AShat(ii,:)'*AShat(ii,:); % (s x s)
% end
% 
% svd_opt = 'svd'; %more accurate, but slower for larger dimensions
% %svd_opt = 'power'; %approximate, faster for larger dimensions
% 
% switch svd_opt
%     case 'svd'
%         [u,sigma,v] = svd(MShat);
%         v1 = u(:,1); %top singular vector of MShat, normalized - s x 1
%     case 'power'
%         v1 = svd_power(MShat);
% end
% 
% v = zeros(n,1);
% v(Shat,1) = v1;
% x_init = phi*v; %ensures that the energy/norm of the initial estimate is close to actual
% x = x_init;
% 
% %initialize x with noisy version of z
x = z+noise;
% %Ground Truth y_p
y_p = (-sign(A*z)+1)/2;

ps = m*0.1;
%p = (-sign(A*x)+1)/2;
%% start descent
fprintf('\n#iter\t\t|y-Ax|\t\t|x-z|\trecovery_prob\n')
for t=1:max_iter 
    
    p = (-sign(A*x)+1)/2;
    %disp('p error before')
     p_err_bef = p-y_p;
%     disp(nnz(p_err_bef))
    %x = cosamp((y_mod-R*p)/sqrt(m), A/sqrt(m), s,100,x); %Its = 10
     %ps = nnz(p_err_bef)+2;
    [x,delta_p] = mod_cosamp(y_mod,p,A,x,R,s,ps);
    
%     indx = find(delta_p);
%     p(indx) = 1-p(indx);
    
    %disp('p_error after')
    %p_err_aft = p-delta_p - y_p;
    %disp(nnz(p_err_aft))
    %err_hist(t+1,1) = norm(y_mod-mod(A*x,R))/norm(y_mod);
    err_hist(t+1,1) = norm((y_mod-y_p*R)-(A*x))/norm(y_mod-y_p*R);
    err_hist(t+1,2) = norm(x-z)/norm(z);
    recovery_prob = nnz(~(p-y_p));
    fprintf('\n%d\t\t%2.8f\t\t%2.4f\t\t%2f\n',t,err_hist(t+1,1),err_hist(t+1,2),recovery_prob)
    if (err_hist(t+1,1) < tol1) | (abs(err_hist(t,2)-err_hist(t+1,2))<tol2)
        break;
    end

end


end