% Testing the robust-ness of the cosamp

clear all;
%Initialize hyper-parameters

n = 1000; %length of the input signal
m = 500; %number of measurements
s = 10; % sparsity
b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
R = 4; %period of the modulo function
iter = 30; %maximum iterations for AltMin based algorithms
iter_gd = 30; %maximum iterations for WF/AF based algorithms (typically take higher no. to convg)
tol1 = 1e-3; %error tolerance for measurements
tol2 = 1e-6;

%Generate the ground truth signal

[z,z_ind] =  generate_signal(n,s,b,1);

%Generate the measurements: y=mod(Ax,R)
[y_mod, y_p, A] = modulo_measure_signal(m,z,R);


%Test the robust-ness of cosamp for perturbed initialization
sig_span = 0.01:0.02:0.2;
err_sig = 0.0*sig_span;
err_init = 0.0*sig_span;

for k = 1:length(sig_span)
    noise = sig_span(k)*randn(n,1);
    [y_mod, y_p, A] = modulo_measure_signal(m,z,R);
    x_init = z + noise;
    
    x = cosamp((y_mod-R*y_p)/sqrt(m), A/sqrt(m),s,10,x_init);
    
    
    %err_init(k) = norm(z-x_init)/norm(z);
    err_sig(k) = norm(x-z)/norm(z)
end

figure;
%plot(sig_span, err_init);
plot(sig_span, err_sig);
title('Mod: Variation of relative error in y with perturbed initialization');

%Test the robust-ness of cosamp for perterbed y_p

perterb_prob = 1:3:30;
err_y_per = 0.0*perterb_prob;
err_x_per = 0.0*perterb_prob;
for i = 1:length(perterb_prob)
    perterb_idx = randperm(s,perterb_prob(i));
    y_perter = y_p;
    y_perter(perterb_idx)= 1 -y_perter(perterb_idx);
    
    x = cosamp((y_mod-R*y_perter)/sqrt(m), A/sqrt(m),s,10);
    
    err_y_per(i) = norm((y_mod-y_p*R)-(y_mod-y_perter*R))/norm((y_mod-y_p*R));
    err_x_per(i) = norm(x-z)/norm(z);
    
end
figure;
plot(perterb_prob,err_y_per);
title('Mod: Variation of relative reconstruction error in x with number of perturbations in p');
figure;
plot(perterb_prob, err_x_per);
title('Mod: Variation of relative error in y with number of perturbations in p');
    
    