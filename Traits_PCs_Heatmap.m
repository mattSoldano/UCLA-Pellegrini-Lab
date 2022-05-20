clear; clc; close;

%% Import traits from metadata

T = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');
traits = T.Properties.VariableNames(3:end);
Tvalues = table2array(T(:,3:end));
Tvalues = knnimpute(Tvalues);

%% Import PCs
PCtable = readtable('pcaData2Comp_ALL_PCs.csv', 'VariableNamingRule', 'preserve');
PCnames = PCtable.Properties.VariableNames;
PCs = table2array(PCtable);


%% Declare X and Y for lasso regression
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

modelcorr = zeros(width(traits), width(PCnames));

for pc = 1:width(PCs)
    for t = 1:width(Tvalues)
        
        Y = Tvalues(:,t);
        X = PCs(:, pc);
        
        %% Determining Correlation

        modelcorr(t, pc) = corr(X,Y);

    end
end

myheatmap = figure(1);

heatmapPCs = heatmap(PCnames, traits, modelcorr);
xlabel('Principal Components');
ylabel('Traits');

saveas(myheatmap, 'PC_Heatmap_noCorrection.png')