% Testing Expectations formulae

%% Testing expectation of sign function: E(a_i*sign(a_i.'*x))
n = 1000; %length of the input signal
%m = 500; %number of measurements
s = 20; % sparsity
b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
R = 4; %period of the modulo function

[z,z_ind] =  generate_signal(n,s,b);

E = zeros(n,1);
n_trials = 1000000.0;

for i=1:n_trials
    a_i = randn(n,1);
    
    E= E + a_i*(sign(a_i.'*z));   
end

E = E/n_trials;

% disp('Expected value is: ')
% disp(E)

disp('difference in the norm is: ')
disp(norm(E-z*sqrt(2.0/pi)))%/norm(z*sqrt(2.0/pi)))

%% Testing first order estimator: E(y_mod*a_i)
n = 1000; %length of the input signal
%m = 500; %number of measurements
s = 20; % sparsity
b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
R = 4; %period of the modulo function

[z,z_ind] =  generate_signal(n,s,b);

E1 = zeros(n,1);
E2 = zeros(n,1);
E3 = zeros(n,1);
E = zeros(n,1);
trialspan = logspace(2,6,5);
err1 = ones(length(trialspan),1);
err2 = err1;
err3 = err1;
err = err1;
errt = err1;
for k = 1:length(trialspan)
    n_trials = trialspan(k)
    for i=1:n_trials
        a_i = randn(n,1); 
        E1 = E1 + a_i*((a_i.'*z));
        E2 = E2 + a_i*(R/2.0);
        E3 = E3 - a_i*(R/2.0)*(sign(a_i.'*z));
        
        E= E + a_i*((a_i.'*z) + (R/2.0)*(1-sign(a_i.'*z)));
    end
    E = E/n_trials;
    E1 = E1/n_trials;
    E2 = E2/n_trials;
    E3 = E3/n_trials;
    
    err1(k) = norm(E1-z)/norm(z);
    err2(k) = norm(E2)/norm(z);
    err3(k) = norm(E3 +(R/2.0)*z*sqrt(2.0/pi))/norm(-(R/2.0)*z*sqrt(2.0/pi));
    
    err(k)=norm(E-(z-(R/2.0)*z*sqrt(2.0/pi)))/norm(z-(R/2.0)*z*sqrt(2.0/pi));
    
    errt(k) = norm(E-(E1+E2+E3));
end
plot(err)
% hold on;
% plot(err1)
% hold on;
% plot(err2)
% hold on;
% plot(err3)


% disp('Expected value is: ')
% disp(E)

% %% Testing second order estimator: E(y^2)
% n = 1000; %length of the input signal
% %m = 500; %number of measurements
% s = 10; % sparsity
% b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
% R = 4; %period of the modulo function
% 
% [z,z_ind] =  generate_signal(n,s,b);
% 
% E = zeros(n,1);
% trialspan = logspace(2,6,5);
% err = ones(length(trialspan),1);
% for k = 1:length(trialspan)
%     n_trials = trialspan(k)
%     for i=1:n_trials
%         a_i = randn(n,1); 
%         E= E + ((a_i.'*z) + (R/2.0)*(1-sign(a_i.'*z))^2)*(a_i.^2);   
%     end
%     E = E/n_trials;
%     [Mg MgS] = sort(E,'descend');
%     
%     S0 = MgS(1:s);
%     Shat = sort(S0); %store indices in sorted order
% 
%     err(k) = sum(ismember(Shat, z_ind'));
%     %err(k)=norm(E-(norm(z)+2*R*R-2*R*((norm(z)^2)/pi)))/norm((norm(z)+2*R*R-2*R*((norm(z)^2)/pi)))
% end
% plot(err)
% 
% %% Testing the first order estimator on real data: E(sum(y_i*a_i))
% clear all;
% n = 1000; %length of the input signal
% mspan = logspace(2,4,3);% 500:1000:10000; %number of measurements
% s = 20; % sparsity
% b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
% R = 4; %period of the modulo function
% 
% [z,z_ind] =  generate_signal(n,s,b);
% 
% 
% err1 = ones(length(mspan),1);
% err2 = ones(length(mspan),1);
% err3 = ones(length(mspan),1);
% norm_x_0 = ones(length(mspan),1);
% for j = 1:length(mspan)
%     m = mspan(j);
%     [y_mod, y_p, A] = modulo_measure_signal(m,z,R);
%     E = zeros(n,1);
%     for i = 1:m
%         E = E + y_mod(i)*A(i,:)';
%     end
%     E = E/m;
%     
%     x_0 = E/ (1-(R/2.0)*(sqrt(2.0/pi)));
%     norm_x_0(j) = norm(x_0);
%     x_0 = x_0/norm(x_0);
%     err3(j)= norm(z-x_0)/norm(z);
%     
%     
%     err1(j) = norm(E - z*(1-(R/2.0)*(sqrt(2.0/pi))))/norm(z*(1-(R/2.0)*(sqrt(2.0/pi))));
%     
%     E = E/norm(E);
%     err2(j)= norm(E-z)/norm(z);
%     
%     
%     
%     
% end
% plot(mspan,err3);
% % plot(err1);
% % hold on;
% % plot(err2);