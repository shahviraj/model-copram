%% Variation of rho (range of A*x) with 'm'
n = 1000;
s = 1000;
b = 1;
mspan= 100:100:500;
n_trials=50;
[z,z_ind] =  generate_signal(n,s,b,1);
rho= zeros(length(mspan),n_trials);
normy= zeros(length(mspan),n_trials);
for j=1:n_trials
        
    for i = 1:length(mspan)

        m = mspan(i);
        A = randn(m,n);
        y = A*z;
        rho(i,j) = max(abs(y));
        %normy(i,j) = norm(y);

    end
  
end
figure;
plot(mspan,mean(rho,2));
%figure;
%plot(mspan,mean(normy,2));