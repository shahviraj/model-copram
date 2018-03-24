clear all;
clc;
disp('starts here:')
%Initialize hyper-parameters
rng ('shuffle')
n = 1000; %length of the input signal
%mspan = 500:500:10000; %number of measurements
%mspan = logspace(2,5,4);

mspan=100:100:1000;
s = 10; % sparsity
b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
R = 1; %period of the modulo function
iter = 5; %maximum iterations for AltMin based algorithms
iter_gd = 30; %maximum iterations for WF/AF based algorithms (typically take higher no. to convg)
tol1 = 1e-4; %error tolerance for measurements
tol2 = 1e-6;
err = zeros(length(mspan),1);
del=1;
max_iter = 30;

supp_recvr=zeros(length(mspan),1);
init_err=zeros(length(mspan),1);
reconst_err=zeros(length(mspan),1);

for j = 1:length(mspan)
    m = mspan(j);
    %Generate the ground truth signal
    [z,z_ind] =  generate_signal(n,s,b);
    
    %z = 10*z;
    %Generate the measurements: y=mod(Ax,R)
    [y_mod, y_p, A] = modulo_measure_signal(m,z,R);
    Marg = zeros(1,n); %marginals
    MShat = zeros(del*s,1); %truncated correlation matrix
    AShat = zeros(m,del*s); %truncated sensing matrix
    supp = zeros(1,n); %indicator for initial support Shat
    y_mod2 = y_mod.^2; %quadratic measurements

    Marg = ((A'.^2)*(y_mod2))/m; % n x 1
    [Mg MgS] = sort(Marg,'descend');
    S0 = MgS(1:del*s); %pick top s-marginals
    Shat = sort(S0); %store indices in sorted order
    
    supp_recvr(j) = sum(ismember(Shat, z_ind'));
    AShat = A(:,Shat);

    card_Marg = m;
    Io = 1:card_Marg;
    
    
    for i = 1:card_Marg
    ii = Io(i);
    MShat = MShat + (y_mod(ii))*AShat(ii,:)'; % (s x s)
    end
    MShat = MShat/card_Marg;
    
    %MShat = MShat/(1-(R/2.0)*sqrt(2.0/pi));
    x_0 = zeros(n,1);
    x_0(Shat,1)= MShat;
    x_0 = x_0/(1-(R/2.0)*sqrt(2.0/pi));
    
    x_0 = norm(z)*(x_0/norm(x_0));
    init_err(j) = norm(z-x_0)/norm(z);
    x= x_0;
    ps = 0.0*m;
    
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
    reconst_err(j) = norm(x-z)/norm(z);
    
    
end
figure;
plot(mspan,reconst_err)
figure;
plot(mspan,init_err)
figure;
plot(mspan,supp_recvr)
%% Verifying for standard compressive sensing

m = 1000;
n = 256;
A =  1/sqrt(m)*randn([m, n]);
x = zeros(n, 1);
k = 8;
supp = randperm(n, k);
x(supp) = randsample([1, -1], k, true);
y = A*x;
% R = 2
% y = y + (1- sign(y))*R/2;
% hist(y)
% M = (A'.^2)*(y.^2)/m;
% sort(M, 'descend');
% [Ms, inds] = sort(M, 'descend');
y = A*x;
M = (A'.^2)*(y.^2)/m;
[Ms, inds] = sort(M, 'descend');


S0 = inds(1:2*k); %pick top s-marginals
Shat = sort(S0); %store indices in sorted order

err = sum(ismember(Shat, supp'))