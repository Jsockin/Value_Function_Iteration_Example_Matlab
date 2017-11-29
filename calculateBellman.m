function [result] = calculateBellman(alpha,beta,delta,psi,z,A,k,kprime,kVals,zPointer,aPointer,...
                                            totalZ,totalA,piZ,piA,valueFun,l1_guess,l1_bar,imagThreshold)
    
    % Pin down kprime on grid of capital in order to find lambda in [0,1]
    if kprime == kVals(length(kVals))
        lowerBound = length(kVals) - 1;
        upperBound = length(kVals);
    else
        lowerBound = max(find(kVals <= kprime));
        upperBound = min(find(kVals > kprime));
    end
    
    % Linearly interpolate value function based on where kprime falls
    lambda          = (kprime - kVals(lowerBound)) / (kVals(upperBound) - kVals(lowerBound));
    interpolatedVF  = valueFun(lowerBound,:,:) + lambda * (valueFun(upperBound,:,:) - valueFun(lowerBound,:,:));
     
    % Calculate optimal labor given k'
    [l1,l2] = solveL1(alpha,delta,psi,z,A,k,kprime,l1_guess,l1_bar);
    labor   = [l1 l2];
    labor   = complex(real(labor),imag(labor).*~(abs(imag(labor))<imagThreshold));

    % Evaluate Bellman equation
    flow                    = utility(alpha,delta,psi,z,A,k,kprime,labor);
    continuation            = 0;
    for nextZ = 1:totalZ
        for nextA = 1:totalA
            continuation    = continuation + (piZ(zPointer,nextZ) * piA(aPointer,nextA)) * interpolatedVF(:,nextZ,nextA);
        end
    end
    bellman                 = flow + beta * continuation;
    
    % Negate bellman for maximization (minimize negative) & save results
    result = -bellman;
    
end

