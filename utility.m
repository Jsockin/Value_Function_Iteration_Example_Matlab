function [util] = utility(alpha,delta,psi,z,A,k,kprime,labor)
    l1 = labor(:,1);
    l2 = labor(:,2);
    first = exp(z) * k^alpha * l1.^(1-alpha) + (1-delta)*k - kprime;
    second = A*l2;
    third = ((l1+l2).^2) / 2;
    util = first.^psi .* second.^(1-psi) - third;
    util(imag(util)~=0) = -10^8;