function [y_clip,y,A] = clip_measure_signal(m,z,R)
%edited 2/15/2017
n = length(z);
%% signal measurement
A = randn(m,n);
y = A*z; %measurements
%y_mod = mod(y,R);
y_clip = y;
y_clip(y_clip<0) = 0;
y_clip(y_clip>R) = R;
end