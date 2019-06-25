function [y_m1, y_m2,y_mp, A_m] = multishot_frwrd(m, z, R)
    
n = length(z);
%% signal measurement
A_m = randn(m/2,n);
y = A_m*z; %measurements

y_m1 = floor((255.99/R)*double((y*(1/2.0) + R/2.0)))/(255.99/R);
%y_m1 = -R + (8/256)*floor(double((y + R)*(255.99/(2*R))));

%y_mod = mod(y,R);
y_mp = (sign(y)-1)/2; %actual phase

%modified modulo function
y_m2 = floor((255.99/R)*double(y - y_mp*R))/(255.99/R);
end