function [result] = runStochastic(alpha,beta,delta,psi,epsilon,kVals,zVals,aVals,piZ,piA,l1_ss,l1_bar,stochDraws,imagThreshold,fminconOptions)

    % Save number of grid options for each state var
    totalK  = length(kVals);
    totalZ  = length(zVals);
    totalA  = length(aVals);
    kmin    = min(kVals);
    kmax    = max(kVals);
    
    % Initiate grid
    [valueFunction,policyK,policyL1,policyL2,policyC1,policyC2] = initializeGrid(totalK,totalZ,totalA,kmin);   

    % Iterate over each state space combo (k,z,A)
    tic
    error = 1;
    totalIterations_stoch = 1;
    while (error > epsilon)

        % Create vector to save results of each draw into
        stochInput  = [randi([1 totalK],stochDraws,1), randi([1 totalZ],stochDraws,1), randi([1 totalA],stochDraws,1)];
        stochKPrime = zeros(stochDraws,1);
        stochVF     = zeros(stochDraws,1);
        stochL1     = zeros(stochDraws,1);
        stochL2     = zeros(stochDraws,1);
        
        priorValueFunction  = valueFunction;
        tempK               = policyK;
        parfor currentDraw = 1:stochDraws

            % Translate current pointers to actual values
            kPointer = stochInput(currentDraw,1);
            k = kVals(kPointer);
            zPointer = stochInput(currentDraw,2);
            z = zVals(zPointer);
            aPointer = stochInput(currentDraw,3);
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

            % Save results for optimal chocie
            stochVF(currentDraw) = -maxVal;
            stochKPrime(currentDraw) = kprimestar;
            stochL1(currentDraw) = l1star;
            stochL2(currentDraw) = l2star;

        end         

        % Save results from stochastic simulations
        for row = 1:stochDraws
            valueFunction(stochInput(row,1),stochInput(row,2),stochInput(row,3)) = stochVF(row);
            policyK(stochInput(row,1),stochInput(row,2),stochInput(row,3))       = stochKPrime(row);
            policyL1(stochInput(row,1),stochInput(row,2),stochInput(row,3))      = stochL1(row);
            policyL2(stochInput(row,1),stochInput(row,2),stochInput(row,3))      = stochL2(row);
        end

        % Calculate new error value and update prior value function estimate
        valueDiff               = valueFunction - priorValueFunction;
        error                   = max(max(max(abs(valueDiff))))
        time_stoch              = toc
        totalIterations_stoch   = totalIterations_stoch + 1;

    end

    % Calculate consumption choices
    [c1_stoch,c2_stoch] = calculateConsumption(alpha,delta,kVals,zVals,aVals,policyK,policyL1,policyL2,policyC1,policyC2);

    % Save results for plotting
    valueF_stoch        = valueFunction;
    k_stoch             = policyK;
    l1_stoch            = policyL1;
    l2_stoch            = policyL2;
    result              = [valueF_stoch k_stoch l1_stoch l2_stoch c1_stoch c2_stoch];
    save('results_stochastic.mat','valueF_stoch','k_stoch','l1_stoch','l2_stoch','c1_stoch','c2_stoch','totalIterations_stoch','time_stoch');

end

