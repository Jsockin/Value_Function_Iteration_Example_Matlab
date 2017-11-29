function [result] = runAccelerator(alpha,beta,delta,psi,epsilon,kVals,zVals,aVals,piZ,piA,l1_ss,l1_bar,imagThreshold,fminconOptions)
    
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
    totalIterations_acc = 1;
    error = 1;
    while (error > epsilon)
        priorValueFunction = valueFunction;
        for kPointer = 1:totalK
            for zPointer = 1:totalZ
                for aPointer = 1:totalA

                    % Translate current pointers to actual values
                    k = kVals(kPointer);
                    z = zVals(zPointer);
                    A = aVals(aPointer);

                    % Maximize only every tenth iteration
                    if mod(totalIterations_acc,10) == 1

                        % Use monotonicity
                        if kPointer > 1
                            kprime0 = policyK(kPointer-1,zPointer,aPointer);
                        else
                            kprime0 = policyK(kPointer,zPointer,aPointer);
                        end

                        % Solve for the optimal kprime
                        bellmanHandle = @(kprime) calculateBellman(alpha,beta,delta,psi,z,A,k,kprime,kVals,zPointer,aPointer,...
                                        totalZ,totalA,piZ,piA,priorValueFunction,l1_ss,l1_bar,imagThreshold);
                        [kprimestar,maxVal] = fmincon(bellmanHandle,kprime0,[],[],[],[],kmin,kmax,[],fminconOptions);

                        % Solve for optimal labor given kprime
                        [l1star,l2star] = solveL1(alpha,delta,psi,z,A,k,kprimestar,l1_ss,l1_bar);

                        % Save results for updating value/policy function
                        currOptimal     = [-maxVal kprimestar l1star l2star];

                    else

                        % Roll forward policy decisions from last iteration
                        accL1       = policyL1(kPointer,zPointer,aPointer);
                        accL2       = policyL2(kPointer,zPointer,aPointer);
                        accKPrime   = policyK(kPointer,zPointer,aPointer);

                        % Pin down accKPrime on grid of capital in order to find lambda in [0,1] 
                        if accKPrime == kVals(length(kVals))
                            lowerBound = length(kVals) - 1;
                            upperBound = length(kVals);
                        else
                            lowerBound = max(find(kVals <= accKPrime));
                            upperBound = min(find(kVals > accKPrime));
                        end

                        % Linearly interpolate value function based on where kprime falls
                        lambda          = (accKPrime - kVals(lowerBound)) / (kVals(upperBound) - kVals(lowerBound));
                        interpolatedVF  = priorValueFunction(lowerBound,:,:) + lambda * (priorValueFunction(upperBound,:,:) - priorValueFunction(lowerBound,:,:));

                        % Calculate updated value function 
                        flow = utility(alpha,delta,psi,z,A,k,accKPrime,[accL1 accL2]);
                        continuation = 0;
                        for nextZ = 1:totalZ
                            for nextA = 1:totalA
                                continuation = continuation + (piZ(zPointer,nextZ) * piA(aPointer,nextA)) * interpolatedVF(:,nextZ,nextA);
                            end
                        end
                        bellman     = flow + beta * continuation;
                        currOptimal = [bellman accKPrime accL1 accL2];

                    end

                    % Save results for optimal choice
                    valueFunction(kPointer,zPointer,aPointer) = currOptimal(1);
                    policyK(kPointer,zPointer,aPointer)       = currOptimal(2);
                    policyL1(kPointer,zPointer,aPointer)      = currOptimal(3);
                    policyL2(kPointer,zPointer,aPointer)      = currOptimal(4);

                end         
            end
        end

        % Calculate new error value and update prior value function estimate
        valueDiff               = valueFunction - priorValueFunction;
        error                   = max(max(max(abs(valueDiff))))
        time_acc                = toc
        totalIterations_acc     = totalIterations_acc + 1;
    end

    % Calculate consumption choices
    [c1_acc,c2_acc] = calculateConsumption(alpha,delta,kVals,zVals,aVals,policyK,policyL1,policyL2,policyC1,policyC2);

    % Save results for testing
    valueF_acc      = valueFunction;
    k_acc           = policyK;
    l1_acc          = policyL1;
    l2_acc          = policyL2;
    result          = [valueF_acc k_acc l1_acc l2_acc c1_acc c2_acc];
    save('results_accelerator.mat','valueF_acc','k_acc','l1_acc','l2_acc','c1_acc','c2_acc','totalIterations_acc','time_acc');

end

