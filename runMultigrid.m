function [result] = runMultigrid(alpha,beta,delta,psi,epsilon,kVals,zVals,aVals,piZ,piA,l1_ss,l1_bar,gridPoints,imagThreshold,fminconOptions)

    % Save number of grid options for each state var
    totalK  = length(kPoints);
    totalZ  = length(zVals);
    totalA  = length(aVals);
    kmin    = min(kVals);
    kmax    = max(kVals);
    
    % Initiate grid
    [valueFunction,policyK,policyL1,policyL2,policyC1,policyC2] = initializeGrid(totalK,totalZ,totalA,kmin);               

    % Run VF iteration with multigrid, linear interpolate between grid
    tic
    error = 1;
    totalIterations_multi = 1;
    while (error > epsilon)
        priorValueFunction  = valueFunction;
        parfor kPointer = 1:totalK
            for zPointer = 1:totalZ
                for aPointer = 1:totalA

                    % Translate current pointers to actual values
                    k = kVals(kPointer);
                    z = zVals(zPointer);
                    A = aVals(aPointer);

                    % Setup for multigrid
                    kprime      = policyK(kPointer,zPointer,aPointer);
                    valueFun    = priorValueFunction;
                    vectorK     = kVals;

                    % Run multigrid
                    for gridIter = 1:length(gridPoints)

                        % Setup parameters based on current grid level
                        if gridIter < length(gridPoints)
                            nextPoints = gridPoints(gridIter+1);
                        else
                            nextPoints = 0;
                        end

                        % Run multigrid algorithm
                        [kprime,vectorK,lowerBound,upperBound,valueFun,maxVal,l1star,l2star] = multigrid(alpha,beta,delta,psi,...
                            z,zPointer,piZ,A,aPointer,piA,k,kprime,nextPoints,vectorK,valueFun,l1_ss,l1_bar,imagThreshold,fminconOptions);    
                    end

                    % Save results for optimal choice
                    valueFunction(kPointer,zPointer,aPointer) = -maxVal;
                    policyK(kPointer,zPointer,aPointer)       = kprime;
                    policyL1(kPointer,zPointer,aPointer)      = l1star;
                    policyL2(kPointer,zPointer,aPointer)      = l2star;

                end
            end
        end

        % Calculate new error value 
        valueDiff               = valueFunction - priorValueFunction;
        error                   = max(max(max(abs(valueDiff))))
        time_multi              = toc
        totalIterations_multi   = totalIterations_multi + 1;
    end

    % Calculate consumption choices
    [c1_multi,c2_multi] = calculateConsumption(alpha,delta,kVals,zVals,aVals,policyK,policyL1,policyL2,policyC1,policyC2);

    % Save results for testing
    valueF_multi        = valueFunction;
    k_multi             = policyK;
    l1_multi            = policyL1;
    l2_multi            = policyL2;
    result              = [valueF_multi k_multi l1_multi l2_multi c1_multi c2_multi];
    save('results_multigrid.mat','valueF_multi','k_multi','l1_multi','l2_multi','c1_multi','c2_multi','totalIterations_multi','time_multi');

end

