function [y_m1, y_m2,y_mp, A_m] = multishot_frwrd(m, z, R, noise)
q = 4096.0;   
n = length(z);
%% signal measurement
A_m = randn(m/2,n);
y = A_m*z; %measurements

y_m1 = floor((q/R)*double((y*(1/2.0) + R/2.0)))/(q/R);

y_m1 = y_m1 + noise(1:m/2,1);
%y_m1 = -R + (8/256)*floor(double((y + R)*(255.99/(2*R))));

%y_mod = mod(y,R);
y_mp = (sign(y)-1)/2; %actual phase

%modified modulo function
y_m2 = floor((q/R)*double(y - y_mp*R))/(q/R);
y_m2 = y_m2 + noise(m/2+1:m,1);
end