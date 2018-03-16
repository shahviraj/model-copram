clear all;
%Initialize hyper-parameters

n = 10000; %length of the input signal
%m = 500; %number of measurements
s = 20; % sparsity
b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
R = 4; %period of the modulo function
iter = 5; %maximum iterations for AltMin based algorithms
iter_gd = 30; %maximum iterations for WF/AF based algorithms (typically take higher no. to convg)
tol1 = 1e-3; %error tolerance for measurements
tol2 = 1e-5;
num_trials=10;
sig_span = 0.0:0.0005:0.005;
m_span = 100:100:500;
%sig_span = 0.04;
err_x = ones(length(m_span), length(sig_span),num_trials, 2);
for j = 1:length(m_span)
    m = m_span(j);
    %Generate the ground truth signal
    [z,z_ind] =  generate_signal(n,s,b);
    
    %Generate the measurements: y=mod(Ax,R)
    [y_mod, y_p, A] = modulo_measure_signal(m,z,R);
    for l = 1:num_trials
        for k = 1:length(sig_span)
            noise = sig_span(k)*randn(n,1);
            err_x(j,k,l,1) = norm(noise)/norm(z);
        % %Use mod_CoPRAM for reconstruction
            [x1,err_hist1,C1] = mod_CoPRAM(y_mod,A,s,iter,R,tol1,tol2,z,noise);

            err_x(j,k,l,2) = norm(x1-z)/norm(z);
        end
    end
end
%Check errors in measurement of p
% recovery_prob = nnz(~(C1-y_p));
% disp('recovered the p correctly for, ')
% disp(recovery_prob)
sum_err_x = sum(err_x,3)/num_trials;
figure;
plot(sum_err_x(1,:,1),sum_err_x(1,:,2));
%legend('m=100');
hold on;
plot(sum_err_x(2,:,1),sum_err_x(2,:,2));
%legend('m=200');
hold on;
plot(sum_err_x(3,:,1),sum_err_x(3,:,2));
%legend('m=300');
hold on;
plot(sum_err_x(4,:,1),sum_err_x(4,:,2));
%legend('m=400');
hold on;
plot(sum_err_x(5,:,1),sum_err_x(5,:,2));
%legend('m=500');
legend('m=100','m=200','m=300','m=400','m=500');
title('Mod: Variation of rel. reconstruction err in x with std. dev. of noise added for initialization, m=500');