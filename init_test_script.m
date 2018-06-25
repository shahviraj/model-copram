clear all;
clc;
%rng ('shuffle')
pr = struct;
%Fixed parameters
pr.n = 1000; %length of the input signal
pr.b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
pr.R = 4; %period of the modulo function
pr.del = 1; %truncation factor for supp estimation
pr.l_pos = 1;
pr.l_neg = -3;
%Tuned parameters
pr.mspan=100:100:1000;
%pr.mspan=8000:8000:8000;
pr.s_span = 3:3:3; % sparsity
pr.amp =1; %amplification factor 
pr.method = 'cosamp';
pr.init_method = 'rcm';
pr.svd_opt = 'svd';

err = zeros(length(pr.mspan),length(pr.s_span));
supp_recvr=zeros(length(pr.mspan),length(pr.s_span));
init_err=zeros(length(pr.mspan),length(pr.s_span));
init_err_mo=zeros(length(pr.mspan),length(pr.s_span));

reconst_err=zeros(length(pr.mspan),length(pr.s_span));

for j = 1:length(pr.mspan)
    m = pr.mspan(j);
    
    for k = 1:length(pr.s_span)
        s = pr.s_span(k);
        %Generate the ground truth signal
        [z,z_ind] =  generate_signal(pr.n,s,pr.b, pr.amp);

%         %Generate the measurements: y=mod(Ax,R)
         [y_mod, y_p, A] = modulo_measure_signal(m,z,pr.R);
        
        %Trying with Phase retrieval measurements
        %[y_abs,y_ph,A] = measure_signal(m,z);
        
        switch pr.init_method
            case 'copram'
                x_0 = copram_init(y_mod,A,s);
            case 'moram'
                x_0 = initial_estimate(A,y_mod,s,pr.R,pr.del,pr.amp);
            case 'raf'
                x_0 = raf_init(A,y_mod,s,pr);
            case 'rcm' %Re-Calculated Measurements
                x_0_mo = initial_estimate(A,y_mod,s,pr.R,pr.del,pr.amp);
                x_0 = rcm_init(A,y_mod,s,pr);
                
        end
        
        %relative error in initial estimate
        init_err(j,k) = norm(z-x_0)/norm(z);
        init_err_mo(j,k) = norm(z-x_0_mo)/norm(z);
    end
    
end


% plot(pr.mspan,init_err(:,1),'-b');
% hold on;
% plot(pr.mspan,init_err_mo(:,1),'-r');

legend_cell = cell(2,1);
legend_cell{1} = 'RCM';
legend_cell{2} = 'MoRAM';
X = [pr.mspan;pr.mspan];
Y = [init_err(:,1)';init_err_mo(:,1)'];
general_plot(X,Y,2,legend_cell,'\textbf{Number of measurements} $\mathbf{(m)}$', '\textbf{Reconstruction Error}; $\mathbf{\frac{||x^*-x||}{||x^*||}}$','\textbf{Comparision of Performance of Initialization Algorithms}')