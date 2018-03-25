clear all;
clc;
rng ('shuffle')

%Fixed parameters
n = 1000; %length of the input signal
b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
tol1 = 1e-5; %error tolerance for measurements
tol2 = 1e-7;
max_iter = 30;
del = 1.0; %amplification factor

%Tuned parameters
mspan=1000:1000:10000;
s_span = 5:5:25; % sparsity
R = 1; %period of the modulo function
del_p=0.05; % ps = del*m (sparsity pertaining to error in p)
method = 'robust-cosamp';


err = zeros(length(mspan),length(s_span));
supp_recvr=zeros(length(mspan),length(s_span));
init_err=zeros(length(mspan),length(s_span));
reconst_err=zeros(length(mspan),length(s_span));

for j = 1:length(mspan)
    m = mspan(j);
    
    for k = 1:length(s_span)
        s = s_span(k);
        %Generate the ground truth signal
        [z,z_ind] =  generate_signal(n,s,b);

        %Generate the measurements: y=mod(Ax,R)
        [y_mod, y_p, A] = modulo_measure_signal(m,z,R);

        %Calculate the initial estimate
        x_0 = initial_estimate(A,y_mod,s,R,del,1);

        init_err(j,k) = norm(z-x_0)/norm(z);

        x= x_0;
        ps = del_p*m;

        fprintf('\n#iter\t\t|y-Ax|\t\t|x-z|\trecovery_prob\n')
        for t=1:max_iter 

            p = (-sign(A*x)+1)/2;

            switch method
                case 'cosamp'
                    x = cosamp((y_mod-R*p)/sqrt(m), A/sqrt(m), s,100,x); %Its = 100
                case 'robust-cosamp'
                    [x,delta_p] = mod_cosamp(y_mod,p,A,x,R,s,ps);
            end

            %err_hist(t+1,1) = norm(y_mod-mod(A*x,R))/norm(y_mod);
            err_hist(t+1,1) = norm((y_mod-y_p*R)-(A*x))/norm(y_mod-y_p*R);
            err_hist(t+1,2) = norm(x-z)/norm(z);
            recovery_prob = nnz(~(p-y_p));
            fprintf('\n%d\t\t%2.8f\t\t%2.4f\t\t%2f\n',t,err_hist(t+1,1),err_hist(t+1,2),recovery_prob)
            if (err_hist(t+1,1) < tol1) | (abs(err_hist(t,2)-err_hist(t+1,2))<tol2)
                break;
            end
        end
    reconst_err(j,k) = norm(x-z)/norm(z);
    end
    
end


construct_plot(reconst_err, mspan,s_span, ['rconst_','r_',num2str(R),'_s_',...
    num2str(s_span(1)),'_',num2str(s_span(end)),'_m_',num2str(mspan(1)),...
    '_',num2str(mspan(end)),'_',method],method);
construct_plot(init_err, mspan,s_span,['init_','r_',num2str(R),'_s_',...
    num2str(s_span(1)),'_',num2str(s_span(end)),'_m_',num2str(mspan(1)),...
    '_',num2str(mspan(end)),'_',method],method);
%figure;
%plot(mspan,supp_recvr)


% %% Verifying for standard compressive sensing
% 
% m = 1000;
% n = 256;
% A =  1/sqrt(m)*randn([m, n]);
% x = zeros(n, 1);
% k = 8;
% supp = randperm(n, k);
% x(supp) = randsample([1, -1], k, true);
% y = A*x;
% % R = 2
% % y = y + (1- sign(y))*R/2;
% % hist(y)
% % M = (A'.^2)*(y.^2)/m;
% % sort(M, 'descend');
% % [Ms, inds] = sort(M, 'descend');
% y = A*x;
% M = (A'.^2)*(y.^2)/m;
% [Ms, inds] = sort(M, 'descend');
% 
% 
% S0 = inds(1:2*k); %pick top s-marginals
% Shat = sort(S0); %store indices in sorted order
% 
% err = sum(ismember(Shat, supp'))