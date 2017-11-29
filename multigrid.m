function [kprime,gridK,lowerBound,upperBound,newVF,maxVal,l1star,l2star] = multigrid(alpha,beta,delta,psi,...
                            z,zPointer,piZ,A,aPointer,piA,k,kprime,nextPoints,vectorK,...
                            valueFun,l1_guess,l1_bar,imagThreshold,fminconOptions)
    
    % Grid setup
    kmin = vectorK(1);
    kmax = vectorK(length(vectorK));
    totalZ = size(valueFun,2);
    totalA = size(valueFun,3);
    
    % Solve for the optimal kprime
    bellmanHandle = @(kprime) calculateBellman(alpha,beta,delta,psi,z,A,k,kprime,vectorK,zPointer,...
                                    aPointer,totalZ,totalA,piZ,piA,valueFun,l1_guess,l1_bar,imagThreshold);
    [kprime,maxVal] = fmincon(bellmanHandle,kprime,[],[],[],[],kmin,kmax,[],fminconOptions);
    
    % Solve for optimal labor given kprime
    [l1star,l2star] = solveL1(alpha,delta,psi,z,A,k,kprime,l1_guess,l1_bar);
    
    % Pin down kprime on grid of capital in order to find lambda in [0,1] 
    if kprime == vectorK(length(vectorK))
        lowerBound = length(vectorK) - 1;
        upperBound = length(vectorK);
    else
        lowerBound = max(find(vectorK <= kprime));
        upperBound = min(find(vectorK > kprime));
    end
    
    % Setup for next layer of grid; otherwise return current values
    if nextPoints > 0
        % Build grid
        gridVals = 1:nextPoints;
        gridIncrement = ((gridVals-1)/(length(gridVals)-1));
        lowerK = vectorK(lowerBound);
        upperK = vectorK(upperBound);
        gridK = (lowerK + gridIncrement*(upperK-lowerK))';
        
        % Update value function for next iteration
        lowerVF = valueFun(lowerBound,:,:);
        upperVF = valueFun(upperBound,:,:);
        newVF = zeros(nextPoints,totalZ,totalA);
        for kNum = 1:nextPoints
            newVF(kNum,:,:) = lowerVF + gridIncrement(kNum)*(upperVF - lowerVF);
        end
    else
        % Return results for last layer of multigrid 
        gridK = vectorK;
        newVF = valueFun;
    end
    
end
