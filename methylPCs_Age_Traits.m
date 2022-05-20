pause
clear; clc; close;
%% Export names
exportName = 'PC_Age_predictTraits.xlsx';

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
Y = traitdata(traitIdx, :);


%% plots

for t = 1:width(methylCompsPCA)
    cornames(t,1) = sprintf("Comp. %d", t);
end

MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
    ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];

[methylCompsVStraits, pval] = corr(Y, scores);

%% Data

for j = 1:width(Vblnames)
mae = zeros(width(scores), 1);
rSq = zeros(width(scores), 1);
agePval = zeros(width(scores), 1);
pcPval = zeros(width(scores), 1);
pc_age_Pval = zeros(width(scores), 1);

for k = 1:width(scores)
    %     X = [Y(:,1), scores(:,k), Y(:,1) .* scores(:,k)];
    X = [Y(:,1), scores(:,k)];
    mdl = fitlm(X , Y(:,j));
    mdlCoefficients = mdl.Coefficients;
    mae(k) = sqrt(mdl.MSE);
    rSq(k) = mdl.Rsquared.Ordinary;
    agePval(k) = mdlCoefficients.pValue(2);
    pcPval(k) = mdlCoefficients.pValue(3);
    %     pc_age_Pval(k) = mdlCoefficients.pValue(4);
end

mdl = fitlm(Y(:,1), Y(:,j));
mdlCoefficients = mdl.Coefficients; 
maeAge = sqrt(mdl.MSE);
rSqAge = mdl.Rsquared.Ordinary;
ageOnlyPval = mdlCoefficients.pValue(2);
%% All Data

ALLmae = [maeAge; mae];
ALLrSq = [rSqAge; rSq];
ALLpcPval = [nan; pcPval];
% ALLpc_age_Pval = [nan; pc_age_Pval];
ALLagePval = [ageOnlyPval; agePval]; 
% allData = [ALLmae, ALLrSq, ALLagePval, ALLpcPval, ALLpc_age_Pval];
allData = [ALLmae, ALLrSq, ALLagePval, ALLpcPval];

%% P-Value, Corr and MAE table
writecell({Vblnames{j}}, exportName, 'Sheet', j, 'Range', 'A1')
writecell({'MAE'}, exportName, 'Sheet', j, 'Range', 'B2')
writecell({'R^2 Values'}, exportName, 'Sheet', j, 'Range', 'C2');
writecell({'Age P Values'}, exportName, 'Sheet', j, 'Range', 'D2')
writecell({'PC P Values'}, exportName, 'Sheet', j, 'Range', 'E2')
% writecell({'PC * Age P Values'}, exportName, 'Sheet', j, 'Range', 'F2')
writecell({'Age Only'}, exportName, 'Sheet', j, 'Range', 'A3');
writecell({cornames{:}}', exportName, 'Sheet', j, 'Range', 'A4')
writematrix(allData, exportName, 'Sheet', j, 'Range', 'B3')
end

%% plots
figs = figure(1);
heatmap(cornames, Vblnames, methylCompsVStraits, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Methyl Principal Component Scores');
ylabel('Traits');

figure(2)
heatmap(cornames, Vblnames, log10(pval), 'ColorLimits', [-4, 1], 'Colormap', MYmap);
xlabel('Methyl Principal Component log10 pValues');
ylabel('Traits');

% subplot(2,2,3)
% plot(mdl)
% xlabel('Age');
% ylabel('Neutrohil %');
% title('fitlm: [Age, PC9], Neutrophil%');