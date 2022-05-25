cd /Users/Sully/Pellegrini/FitnessData
pause
clear; clc; close;

%% bring in data

%%%%%% Need to remove outliers

load('bloodMethylPC_Trait_Associations');
cd BloodPCresults
valueTable=readtable('PC_Trait_AssociationModel_Blood.xlsx', 'VariableNamingRule', 'preserve');
%% plots

for i = 1:10
    trait = valueTable.Var11{i};
    PC = valueTable.Var1(i);
    
    for si = 1:height(predictY)
        if strcmp(predictY{si,1}, trait) && isequal(predictY{si,2}, PC)
            pY = cell2mat(predictY(si,3));
            % predicted trait values
            break
        end
    end
    
    for ai = 1:width(Vblnames)
        if strcmp(Vblnames{ai}, trait)
            aY = Y(:,ai);
            break
        end
    end
    
    
    fig = figure(1);
    
%     aY = aY(l3d.(trait), :);
    % actual trait values --> l3d was to remove outliers

    poly = polyfit(aY, pY, 1);
    fit = polyval(poly, aY);
    plot(aY, pY, 'o', aY, fit, '-');
    title(sprintf('PC %.0f', PC));
    xlabel(trait);
    ylabel(['Predicted ', trait]);
    
    savename = sprintf('%s_bloodPC_%.0f_linearfit.png', trait, PC);
    saveas(fig, savename);
    close all
    clear fig
end