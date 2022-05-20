pause
clear; clc; close;

%% bring in data

load('SiteSpecificAssocationsRemovedOutliers');
valueTable=readtable('SitePredict_removeOutliers.xlsx', 'VariableNamingRule', 'preserve');
%% plots

for i = 1:10
    trait = valueTable.Traits{i};
    site = valueTable.Site{i};
    
    for si = 1:height(predictY)
        if strcmp(predictY{si,1}, trait) && strcmp(predictY{si,2}, site)
            pY = cell2mat(predictY(si,3));
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
    
    scatter(aY, pY);
    title(site);
    xlabel(trait);
    ylabel(['Predicted ', trait]);
    
    savename = [trait, '_', site, '.png'];
    saveas(fig, savename);
    close all
    clear fig
end