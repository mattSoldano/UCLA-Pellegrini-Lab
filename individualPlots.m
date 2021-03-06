pause
clear; clc; close;

%% bring in data

load('bloodMethylPC_Trait_Associations');
valueTable=readtable('PC_Trait_AssociationModel_Blood.xlsx', 'VariableNamingRule', 'preserve');
%% plots

for i = 1:10
    trait = valueTable.Traits{i};
    site = valueTable.Site{i};
    
    for si = 1:height(predictY)
        if strcmp(predictY{si,1}, trait) && strcmp(predictY{si,2}, site)
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
    
    site = regexprep(site, '_', ':');
    
    fig = figure(1);
    
    aY = aY(l3d.(trait), :);
    % actual trait values

    poly = polyfit(aY, pY, 1);
    fit = polyval(poly, aY);
    plot(aY, pY, 'o', aY, fit, '-');
    title(site);
    xlabel(trait);
    ylabel(['Predicted ', trait]);
    
    savename = [trait, '_', site, '.png'];
    saveas(fig, savename);
    close all
    clear fig
end