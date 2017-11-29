function [l1,l2] = solveL1(alpha,delta,psi,z,A,k,kprime,l1,l1_bar)

    % Setup for bisection algorithm to pin down l1 given l2
    lowerL1 = 0;
    upperL1 = l1_bar;
    difference = 1;
    epsilon = 10^-6;
    imagThreshold = 10^-4;
    
    % Solve for optimal l1 given l1 using bisection
    while abs(difference) > epsilon
        c1 = exp(z)*k^alpha*l1^(1-alpha) + (1-delta)*k - kprime;
        l2 = ((1-psi) * c1) / (psi * (1-alpha) * exp(z) * (k/l1)^alpha);
        lhs = (l1 + l2)^(1/psi) * l1^alpha;
        rhs = psi * (A*(1-psi))^((1-psi)/psi) * (1-alpha) * exp(z) * k^alpha;
        difference = lhs - rhs;
        
        % Value for l1 from bisection leads to negative consumption
        if abs(imag(difference)) > imagThreshold
            l1 = -1;
            break;
        elseif difference > 0
            upperL1 = l1;
            l1 = l1 - 0.5 * (l1 - lowerL1);
        elseif difference < 0
            lowerL1 = l1;
            l1 = l1 + 0.5 * (upperL1 - l1);
        end

    end
    
end

