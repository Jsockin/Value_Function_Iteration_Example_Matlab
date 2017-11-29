function [] = plotFunctions(endFileName,endVarName,kVals,zVals,aVals)

    % Load relevant variables
    fileName = ['results_' endFileName '.mat'];
    load(fileName);
    % Set parameters / values for loop plotting
    listResults = {['valueF_' endVarName],['c1_' endVarName],['c2_' endVarName],...
                   ['k_' endVarName],['l1_' endVarName],['l2_' endVarName]};
    listNames = {'v(k,z,' , 'c_1(k,z,' , 'c_2(k,z,', 'kprime(k,z,' , 'l_1(k,z,' , 'l_2(k,z,'};
    totalVars = length(listResults);
    yLimits = zeros(totalVars,2);
    for varNum = 1:totalVars
        yLimits(varNum,1) = 0.95 * min(min(min(eval(listResults{varNum}))));
        yLimits(varNum,2) = 1.05 * max(max(max(eval(listResults{varNum}))));
    end
    legendString = strings(length(zVals),1);
    for zNum = 1:length(zVals)
        legendString(zNum,1) = ['z = ' num2str(zVals(zNum))];
    end

    % Plot the value function iteration results
    rows = 3;
    cols = 3;
    for plotRow = 1:totalVars

        if mod(plotRow,rows) == 1
            figVFIter = figure(plotRow);
        end

        currPlotMat = eval(listResults{plotRow});  

        for plotCol = 1:cols

            if mod(plotRow,3) == 0
                jump = mod(plotRow-1,3);
            else
                jump = mod(plotRow,3)-1;
            end
            subplot(rows,cols,jump*cols+plotCol);
            plot(kVals,currPlotMat(:,:,plotCol));
            ylim([yLimits(plotRow,1) yLimits(plotRow,2)]);
            if plotCol == 1 && mod(plotRow,rows) == 1
                plotLegend = legend(legendString,'Location','Best','FontSize',5);
                set(plotLegend,'Box','off','color','none');
            end
            title([char(listNames(plotRow)) num2str(aVals(plotCol)) ')']);

        end

        if mod(plotRow,rows) == 0
            saveas(figVFIter,['Figures/Plots_vfi_' endVarName '_' num2str(plotRow) '.png']);
            close(figVFIter);
        end

    end

end