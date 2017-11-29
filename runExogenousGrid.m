function [result] = runExogenousGrid(alpha,beta,delta,psi,epsilon,kVals,zVals,aVals,piZ,piA,l1_ss,l1_bar,imagThreshold,fminconOptions)

    % Save number of grid options for each state var
    totalK = length(kPoints);
    totalZ = length(zVals);
    totalA = length(aVals);
    kmin    = min(kVals);
    kmax    = max(kVals);

    % Initiate grid
    [valueFunction,policyK,policyL1,policyL2,policyC1,policyC2] = initializeGrid(totalK,totalZ,totalA,kmin);               

    % Iterate over each state space combo (k,z,A)
    tic
    error = 1;
    totalIterations_exogrid = 1;
    while (error > epsilon)
        priorValueFunction  = valueFunction;
        tempK = policyK;
        parfor kPointer = 1:totalK
            for zPointer = 1:totalZ
                for aPointer = 1:totalA

                    % Translate current pointers to actual values
                    k = kVals(kPointer);
                    z = zVals(zPointer);
                    A = aVals(aPointer);

                    % Use monotonicity
                    if kPointer > 1
                        kprime0 = tempK(kPointer-1,zPointer,aPointer);
                    else
                        kprime0 = tempK(kPointer,zPointer,aPointer);
                    end

                    % Solve for the optimal kprime
                    bellmanHandle = @(kprime) calculateBellman(alpha,beta,delta,psi,z,A,k,kprime,kVals,zPointer,aPointer,...
                                        totalZ,totalA,piZ,piA,priorValueFunction,l1_ss,l1_bar,imagThreshold);
                    [kprimestar,maxVal] = fmincon(bellmanHandle,kprime0,[],[],[],[],kmin,kmax,[],fminconOptions);

                    % Solve for optimal labor given kprime
                    [l1star,l2star] = solveL1(alpha,delta,psi,z,A,k,kprimestar,l1_ss,l1_bar);

                    % Save results for optimal choice
                    valueFunction(kPointer,zPointer,aPointer) = -maxVal;
                    policyK(kPointer,zPointer,aPointer)       = kprimestar;
                    policyL1(kPointer,zPointer,aPointer)      = l1star;
                    policyL2(kPointer,zPointer,aPointer)      = l2star;

                end         
            end
        end

        % Calculate new error value and update prior value function estimate
        valueDiff               = valueFunction - priorValueFunction;
        error                   = max(max(max(abs(valueDiff))))
        time_exogrid            = toc
        totalIterations_exogrid = totalIterations_exogrid + 1;
        
    end

    % Calculate consumption choices
    [c1_exogrid,c2_exogrid] = calculateConsumption(alpha,delta,kVals,zVals,aVals,policyK,policyL1,policyL2,policyC1,policyC2);

    % Save results 
    valueF_exogrid          = valueFunction;
    k_exogrid               = policyK;
    l1_exogrid              = policyL1;
    l2_exogrid              = policyL2;
    result                  = [valueF_exogrid k_exogrid l1_exogrid l2_exogrid c1_exogrid c2_exogrid];
    save('results_exogrid.mat','valueF_exogrid','k_exogrid','l1_exogrid','l2_exogrid','c1_exogrid','c2_exogrid','totalIterations_exogrid','time_exogrid');

end

