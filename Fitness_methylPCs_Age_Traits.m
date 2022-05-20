pause
clear; clc; close;
%% Export names
exportName = 'Summary_Fitness_PC_Age_Gender_predictTraits.xlsx';

%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
methylT = readtable('MethMatrixCG-HumanBloodFitness96manual_20-100.csv', 'VariableNamingRule', 'preserve');
methylData = table2array(methylT(:, 2:end));
methylSites = methylT.Properties.VariableNames(2:end);
samplenames = methylT.sampleId;
%% Import traits and import the sample names we have methylation values for

T = readtable('Traits_ProsperProteome103_120821mjt.csv', 'VariableNamingRule', 'preserve');
Vblnames = T.Properties.VariableNames(2:end);
traitdata = table2array(T(:, 2:end));
traitnames = T.sampleID;

%% Compare names from metadata and methylation value names to make sure they are the same
traitIdx = zeros(length(samplenames),1);

for kit=1:numel(samplenames)
    if length(find(traitnames == samplenames(kit))) == 1
        traitIdx(kit) = find(traitnames == samplenames(kit));
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

age = Y(:,3);
gender = Y(:,2);

for j = 1:width(Vblnames)
    Vblnames{j}
    mae = zeros(width(scores), 1);
    rSq = zeros(width(scores), 1);
    agePval = zeros(width(scores), 1);
    pcPval = zeros(width(scores), 1);
    genderPval = zeros(width(scores), 1);
    pc_age_Pval = zeros(width(scores), 1);
    
    for k = 1:width(scores)
        %     X = [Y(:,1), scores(:,k), Y(:,1) .* scores(:,k)];
        X = [age, gender, scores(:,k)];
        mdl = fitlm(X , Y(:,j));
        mdlCoefficients = mdl.Coefficients;
        mae(k) = sqrt(mdl.MSE);
        rSq(k) = mdl.Rsquared.Ordinary;
        agePval(k) = log10(mdlCoefficients.pValue(2));
        genderPval(k) = log10(mdlCoefficients.pValue(3));
        pcPval(k) = log10(mdlCoefficients.pValue(4));
        %     pc_age_Pval(k) = mdlCoefficients.pValue(4);
    end
    
    mdl = fitlm(age, Y(:,j));
    mdlCoefficients = mdl.Coefficients;
    maeAge = sqrt(mdl.MSE);
    rSqAge = mdl.Rsquared.Ordinary;
    ageOnlyPval = log10(mdlCoefficients.pValue(2));
    
    mdl = fitlm(gender, Y(:,j));
    mdlCoefficients = mdl.Coefficients;
    maeGender = sqrt(mdl.MSE);
    rSqGender = mdl.Rsquared.Ordinary;
    GenderOnlyPval = log10(mdlCoefficients.pValue(2));
    %% All Data
    
    ALLmae = [maeAge; maeGender; mae];
    ALLrSq = [rSqAge; rSqGender; rSq];
    ALLpcPval = [nan; nan; pcPval];
    ALLagePval = [ageOnlyPval; nan; agePval];
    ALLgenderPval = [nan; GenderOnlyPval; genderPval];
    allData = [ALLmae, ALLrSq, ALLagePval, ALLgenderPval, ALLpcPval];
    
    if numel(find(allData == -inf)) >= 1 || numel(find(allData == inf)) >= 1
        allData = num2cell(allData);
        allData = cellfun(@num2str, allData, 'UniformOutput', false);
        writecell(allData, exportName, 'Sheet', j, 'Range', 'B3')
    else
        writematrix(allData, exportName, 'Sheet', j, 'Range', 'B3')
    end
    
    %% P-Value, Corr and MAE table
   
    
end

%% plots
% figs = figure(1);
% heatmap(cornames, Vblnames, methylCompsVStraits, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
% xlabel('Methyl Principal Component Scores');
% ylabel('Traits');
%
% figure(2)
% heatmap(cornames, Vblnames, log10(pval), 'ColorLimits', [-4, 1], 'Colormap', MYmap);
% xlabel('Methyl Principal Component log10 pValues');
% ylabel('Traits');

% subplot(2,2,3)
% plot(mdl)
% xlabel('Age');
% ylabel('Neutrohil %');
% title('fitlm: [Age, PC9], Neutrophil%');