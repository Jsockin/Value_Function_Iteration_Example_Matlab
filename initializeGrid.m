function [vf,k,l1,l2,c1,c2] = initializeGrid(totalK,totalZ,totalA,kmin)
    
    % Initalize policy functions
    k  = kmin*ones(totalK,totalZ,totalA);
    l1 = zeros(totalK,totalZ,totalA);
    l2 = zeros(totalK,totalZ,totalA);
    c1 = zeros(totalK,totalZ,totalA);
    c2 = zeros(totalK,totalZ,totalA);
    
    % Initialize value function with curvature
    vf = zeros(totalK,totalZ,totalA);
    for row = 1:totalK
        vf(row,:,:)  = log(row);
    end
    
end

