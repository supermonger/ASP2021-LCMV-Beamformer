function W = LCMV(R,theta_s,theta_i,gi,N)
    C = zeros(N,2);
    C(:,1) = a_of_theta(theta_s,N);
    C(:,2) = a_of_theta(theta_i,N);
    g = [1;gi];
    W = (R\C)/(C'/R*C)*g;
end