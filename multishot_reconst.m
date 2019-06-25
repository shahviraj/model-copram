function x = multishot_reconst(y_m1, y_m2, A_m, R, opts)
%function [E2,k2] = multishot_reconst(y_m1,y_m2,R)    

    m = 2*length(y_m1);
    
    k1 = zeros(m/2,1);

    E1 = y_m1 + k1*(R);

    k2 = floor(double(2*E1)/(R));

    E2 = y_m2 + k2*(R); 
    
    E2 = E2  - R;
    
    x = spg_bp(A_m/sqrt(m/2.0),E2/sqrt(m/2.0),opts);
    
  