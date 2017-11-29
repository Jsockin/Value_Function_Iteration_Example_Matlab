function [] = plotComparison(endFileName1,endVarName1,endFileName2,endVarName2,kVals1,kVals2,zVals,aVals)

    % Load relevant variables
    fileName1 = ['results_' endFileName1 '.mat'];
    load(fileName1);
    fileName2 = ['results_' endFileName2 '.mat'];
    load(fileName2);
    
    % Fix z and A for comparison (at steady state)
    zPointer = 3;
    aPointer = 2;
    
    % Set parameters / values for loop plotting
    listResults1 = {['valueF_' endVarName1],['c1_' endVarName1],['c2_' endVarName1],...
                   ['k_' endVarName1],['l1_' endVarName1],['l2_' endVarName1]};
    listResults2 = {['valueF_' endVarName2],['c1_' endVarName2],['c2_' endVarName2],...
                   ['k_' endVarName2],['l1_' endVarName2],['l2_' endVarName2]};
    listNames = {'v(k,' , 'c_1(k,' , 'c_2(k,', 'kprime(k,' , 'l_1(k,' , 'l_2(k,'};
    totalVars = length(listResults1);
    legendString = strings(2,1);
    legendString(1,1) = endFileName1;
    legendString(2,1) = endFileName2;
    
    % Plot the value function iteration results
    rows = 3;
    cols = 2;
    figVFIter = figure(1);
    for varNum = 1:totalVars
        
        % Load plot data for comparison
        currMat1 = eval(listResults1{varNum}); 
        plotData1 = currMat1(:,zPointer,aPointer);  
        currMat2 = eval(listResults2{varNum});  
        plotData2 = currMat2(:,zPointer,aPointer);

        subplot(rows,cols,varNum);
        plot(kVals1,plotData1);
        hold on;
        plot(kVals2,plotData2);
        yLowerLim = 0.98 * min(min(min(min(plotData1))),min(min(min(plotData1))));
        yUpperLim = 1.02 * max(max(max(max(plotData2))),max(max(max(plotData2))));
        ylim([yLowerLim yUpperLim]);
        if varNum == 1
            plotLegend = legend(legendString,'Location','Best','FontSize',5);
            set(plotLegend,'Box','off','color','none');
        end
        title([char(listNames(varNum)) num2str(zVals(zPointer)) ',' num2str(aVals(aPointer)) ')']);

    end
    
    saveas(figVFIter,['Figures/Plots_compare_' endVarName1 '_' endVarName2 '.png']);
    close(figVFIter);
    
end