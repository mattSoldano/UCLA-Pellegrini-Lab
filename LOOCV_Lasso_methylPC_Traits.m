pause
clear; clc; close;
%% Export names
exportName = 'PC_Age_predictTraits_lasso.xlsx';

%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
methylData = dlmread('MethMatrixTurk.csv', ',');
methylCoords = readcell('MethMatrixProbeNames.csv');
methylCoords = methylCoords(2:end);
samplenames = readcell('MethMatrixSampleNames.csv');
samplenames = (erase(samplenames(2:end), 'R'));
samplenames = regexprep(samplenames, '-', '_');
%% Import traits and import the sample names we have methylation values for

T = readtable('Turkish_Buccal_Traits2.csv', 'VariableNamingRule', 'preserve');
Vblnames = T.Properties.VariableNames;

for g = 1:height(T.Gender)
    if isequal(T.Gender(g), "M") == true
        T.Gender{g} = 0;
    else
        T.Gender{g} = 1;
    end
end
% change the gender value --> Female = 1, male = 0
T.Gender = cell2mat(T.Gender);
desiredTraits = [3 4 8:29 31 32];
traitnames = readcell('MethMatrixTraitNames.csv');
traitnames = traitnames(2:end);
traitdata = table2array(T(:, desiredTraits));
Vblnames = Vblnames(desiredTraits);

%% Remove sample C40 because of NaN values
C40t = find(strcmp(traitnames, "C40"));
traitnames(C40t) = [];
traitdata(C40t,:) = [];

C40m = find(strcmp(samplenames, "Turk_C40"));
methylData(C40m, :) = [];
samplenames(C40m, :) = [];

%% Rename traitnames
turk = 'Turk_';
for s = 1:height(traitnames)
    traitnames{s} = append(turk, num2str(traitnames{s}));
end

%% Compare names from metadata and methylation value names to make sure they are the same
traitIdx = zeros(length(samplenames),1);

for kit=1:numel(samplenames)
    if length(find(strcmp(traitnames,samplenames{kit}))) == 1
        traitIdx(kit) = find(strcmp(traitnames,samplenames{kit}));
    end
end
%% PCA

[methylCompsPCA, scores] = pca(methylData);
X = [traitdata(traitIdx, 1), scores];
%% Data

uniq = unique(samplenames);
numUniqueSamples = numel(uniq);
numsamples = height(X);
for traitid = 1:width(desiredTraits)
    traitid
    Y = traitdata(traitIdx, traitid);
    trait = Vblnames{traitid};
    trait = regexprep(trait, ' ', '_');
    trait = regexprep(trait, '-', '_');
    trait = regexprep(trait, '%', 'percent');
    trait = regexprep(trait, '#', 'number');
    clear B; clear YTrain; clear XTrain; clear idxLambda1SE; clear FitInfo
    pause(90);
    for k=1:numUniqueSamples
        idx = [1:numsamples]';
        id = table(idx, samplenames);
        clear testSample
        for j = 1:numsamples
            if isequal(uniq{k},id.samplenames{j}) == true
                testSample(1) = j;
                id(j,:) = [];
                for sa = 1:height(id.samplenames)
                    if isequal(uniq{k},id.samplenames{sa})
                        testSample(2) = sa + 1;
                        id(sa,:) = [];
                        break
                    end
                end
                break
            else
            end
        end
        
        XTrain = X(id.idx, :);
        YTrain = Y(id.idx, :);
        [B,FitInfo] = lasso(XTrain,YTrain,'Alpha',1.0,'CV',10);
        
        idxLambda1SE = FitInfo.Index1SE;
        coef(:,k) = B(:,idxLambda1SE);
        coef0(:,k) = FitInfo.Intercept(idxLambda1SE);
        
        for dups = 1:numel(testSample)
            yhat(testSample(dups), traitid) = X(testSample(dups),:)*coef(:,k) + coef0(:,k);
        end
    end
    [modelcorr(traitid), pVals(traitid)] = corr(yhat(:, traitid),Y);
    mae(traitid) = median(abs(Y(:,1) - yhat(:,traitid)));
    coefTraits.(trait) = coef;
    coef0Traits.(trait) = coef0;
end
save('AllPC_TraitLasso');
%% Heatmap

MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
    ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];

Y = traitdata(traitIdx, :);

for t = 1:width(methylCompsPCA)
    cornames(t,1) = sprintf("Comp. %d", t);
end

% Try this again with rescaled data after you output all the crucial info
t_pt = figure(1);
heatmap(Vblnames, Vblnames, corr(yhat, Y), 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Traits');
ylabel('Predicted Traits');

%% Save the data

writecell({'R Values'}, exportName, 'Range', 'B1');
writecell({Vblnames{:}}', exportName, 'Range', 'A2')
writematrix(modelcorr', exportName, 'Range', 'B2')
writecell({'log10(P-Values)'}, exportName, 'Range', 'C1')
writematrix(log10(pVals)', exportName, 'Range', 'C2')
writecell({'MAE:'}, exportName, 'Range', 'D1')
writematrix(mae', exportName, 'Range', 'D2')

j = 1;

for sheet = 2:(width(desiredTraits)+1)
    
    trait = Vblnames{j};
    trait = regexprep(trait, ' ', '_');
    trait = regexprep(trait, '-', '_');
    trait = regexprep(trait, '%', 'percent');
    trait = regexprep(trait, '#', 'number');
    
    writecell({'Coefficient Values'}, exportName, 'Sheet', sheet, 'Range', 'A1');
    writecell({trait}, exportName, 'Sheet', sheet, 'Range', 'C1');
    writecell({'Age'}, exportName, 'Sheet', sheet, 'Range', 'A2')
    writecell({cornames{:}}', exportName, 'Sheet', sheet, 'Range', 'A3')
    writematrix(coefTraits.(trait), exportName, 'Sheet', sheet, 'Range', 'B2')

    j = j + 1;
end