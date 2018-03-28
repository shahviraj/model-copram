clear all;
clc;
rng ('shuffle')
pr = struct;
%Fixed parameters
pr.n = 1000; %length of the input signal
pr.b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
pr.tol1 = 1e-5; %error tolerance for measurements
pr.tol2 = 1e-7;
pr.max_iter = 30;
pr.del = 1.0; %amplification factor

%Tuned parameters
pr.mspan1 = [100:100:1000];
pr.mspan2 = [1000:1000:10000];
pr.mspan=[pr.mspan1,pr.mspan2];
%pr.mspan=5000:1000:8000;
pr.s_span = 5:5:5; % sparsity
pr.R = 1; %period of the modulo function
pr.del_p=0.05; % ps = del*m (sparsity pertaining to error in p)
pr.method = 'cosamp';


err = zeros(length(pr.mspan),length(pr.s_span));
supp_recvr=zeros(length(pr.mspan),length(pr.s_span));
init_err=zeros(length(pr.mspan),length(pr.s_span));
reconst_err=zeros(length(pr.mspan),length(pr.s_span));

for j = 1:length(pr.mspan)
    m = pr.mspan(j);
    
    for k = 1:length(pr.s_span)
        s = pr.s_span(k);
        %Generate the ground truth signal
        [z,z_ind] =  generate_signal(pr.n,s,pr.b);

        %Generate the measurements: y=mod(Ax,R)
        [y_mod, y_p, A] = modulo_measure_signal(m,z,pr.R);

        %Calculate the initial estimate
        x_0 = initial_estimate(A,y_mod,s,pr.R,pr.del,1);
        
        %relative error in initial estimate
        init_err(j,k) = norm(z-x_0)/norm(z);
        
        %Alt-Min
        x= x_0;
        ps = pr.del_p*m;

        fprintf('\n#iter\t\t|y-Ax|\t\t|x-z|\trecovery_prob\n')
        for t=1:pr.max_iter 

            p = (-sign(A*x)+1)/2;

            switch pr.method
                case 'cosamp'
                    x = cosamp((y_mod-pr.R*p)/sqrt(m), A/sqrt(m),s,100,x); %Its = 100
                case 'robust-cosamp'
                    [x,delta_p] = mod_cosamp(y_mod,p,A,x,pr.R,s,ps);
            end

            %err_hist(t+1,1) = norm(y_mod-mod(A*x,R))/norm(y_mod);
            err_hist(t+1,1) = norm((y_mod-y_p*pr.R)-(A*x))/norm(y_mod-y_p*pr.R);
            err_hist(t+1,2) = norm(x-z)/norm(z);
            recovery_prob = nnz(~(p-y_p));
            fprintf('\n%d\t\t%2.8f\t\t%2.4f\t\t%2f\n',t,err_hist(t+1,1),err_hist(t+1,2),recovery_prob)
%             if (err_hist(t+1,1) < pr.tol1) | (abs(err_hist(t,2)-err_hist(t+1,2))<pr.tol2)
%                 break;
%             end
        end
        % Relative reconstruction error
        reconst_err(j,k) = norm(x-z)/norm(z);
    end
    
end


construct_subplots(reconst_err,pr,['rconst_','r_',num2str(pr.R),'_s_',...
    num2str(pr.s_span(1)),'_',num2str(pr.s_span(end)),'_m_',num2str(pr.mspan(1)),...
    '_',num2str(pr.mspan(end)),'_',pr.method],pr.method);
% construct_plot(init_err,pr,['init_','r_',num2str(pr.R),'_s_',...
%     num2str(pr.s_span(1)),'_',num2str(pr.s_span(end)),'_m_',num2str(pr.mspan(1)),...
%     '_',num2str(pr.mspan(end)),'_',pr.method],pr.method);