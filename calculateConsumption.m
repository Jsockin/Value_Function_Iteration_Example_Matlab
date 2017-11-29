function [c1,c2] = calculateConsumption(alpha,delta,kVals,zVals,aVals,policyK,policyL1,policyL2,policyC1,policyC2)

    % Calculate consumption choices
    c1 = policyC1;
    c2 = policyC2;
    for zPointer = 1:length(zVals)
        for aPointer = 1:length(aVals)
            c1(:,zPointer,aPointer) = exp(zVals(zPointer)) .* kVals.^alpha .* policyL1(:,zPointer,aPointer).^(1-alpha) ...
                                                + (1-delta).*kVals - policyK(:,zPointer,aPointer);
            c2(:,zPointer,aPointer) = aVals(aPointer).*policyL2(:,zPointer,aPointer);
        end
    end

end

