pause
clear; clc; close;
%% Export names
exportName = 'Fitness_PC_Age_Gender_predictTraits_LOOCV.xlsx';

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
    cornames(t,1) = sprintf("PC %d", t);
end

MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
    ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];

[methylCompsVStraits, pval] = corr(Y, scores);

%% Data

age = Y(:,3);
gender = Y(:,2);

numsamples = height(Y);
n = 1;

for j = 1:width(Vblnames)
    trait = Vblnames{j}
    trait = regexprep(trait, ' ', '_');
    trait = regexprep(trait, '-', '_');
    trait = regexprep(trait, '%', 'percent');
    trait = regexprep(trait, '#', 'number');
    trait = regexprep(trait, '\.', '_');
    mae = zeros(width(scores), 1);
    rSq = zeros(width(scores), 1);
    agePval = zeros(width(scores), 1);
    pcPval= zeros(width(scores), 1);
    genderPval = zeros(width(scores), 1);
    
    for k = 1:width(scores)
        X = [age, gender, scores(:,k)];
        
        for loo=1:numsamples
            idx = [1:numsamples];
            idx(loo) = [];
            
            XTrain = X(idx, :);
            YTrain = Y(idx, j);
            LOOmdl = fitlm(XTrain , YTrain);
            
            LOOmdlCoefficients = LOOmdl.Coefficients;
            mae(k) = sqrt(LOOmdl.MSE);
            rSq(k) = LOOmdl.Rsquared.Ordinary;
            ageEST(k) = LOOmdlCoefficients.Estimate(2);
            genderEST(k) = LOOmdlCoefficients.Estimate(3);
            pcEST(k) = LOOmdlCoefficients.Estimate(4);
            agePval(k) = log10(LOOmdlCoefficients.pValue(2));
            genderPval(k) = log10(LOOmdlCoefficients.pValue(3));
            pcPval(k) = log10(LOOmdlCoefficients.pValue(4));
        end
        
        averageMAE(k, j) = mean(mae);
        stdMAE(k, j) = std(mae);
        averageRSQ(k, j) = mean(rSq);
        stdRSQ(k, j) = std(rSq);
        averageAgeP(k, j) = mean(agePval);
        stdAgeP(k, j) = std(agePval);
        averageGenderP(k, j) = mean(genderPval);
        stdGenderP(k, j) = std(genderPval);
        averagePCP(k, j) = mean(pcPval);
        stdPCP(k, j) = std(pcPval);
        averageAgeE(k, j) = mean(ageEST);
        stdAgeE(k, j) = std(ageEST);
        averageGenderE(k, j) = mean(genderEST);
        stdGenderE(k, j) = std(genderEST);
        averagePCE(k, j) = mean(pcEST);
        stdPCE(k, j) = std(pcEST);
        
    end
    
    
    
    %     mdl = fitlm(age, yhat(:,j));
    %     mdlCoefficients = mdl.Coefficients;
    %     maeAge = sqrt(mdl.MSE);
    %     rSqAge = mdl.Rsquared.Ordinary;
    %     ageOnlyPval = log10(mdlCoefficients.pValue(2));
    %
    %     mdl = fitlm(gender, yhat(:,j));
    %     mdlCoefficients = mdl.Coefficients;
    %     maeGender = sqrt(mdl.MSE);
    %     rSqGender = mdl.Rsquared.Ordinary;
    %     GenderOnlyPval = log10(mdlCoefficients.pValue(2));
    %% All Data
    
    %     ALLmae = [maeAge; maeGender; mae];
    %     ALLrSq = [rSqAge; rSqGender; rSq];
    %     ALLpcPval = [nan; nan; pcPval];
    %     ALLagePval = [ageOnlyPval; nan; agePval];
    %     ALLgenderPval = [nan; GenderOnlyPval; genderPval];
    %     exportSamplenames = [samplenames; nan];
    %     allData = [ALLmae, ALLrSq, ALLagePval, ALLgenderPval, ALLpcPval, exportSamplenames];
    %
    %     if numel(find(allData == -inf)) >= 1 || numel(find(allData == inf)) >= 1
    %         allData = num2cell(allData);
    %         allData = cellfun(@num2str, allData, 'UniformOutput', false);
    %         writecell(allData, exportName, 'Sheet', j, 'Range', 'B3')
    %     else
    %         writematrix(allData, exportName, 'Sheet', j, 'Range', 'B3')
    %     end
    %
    %     %% P-Value, Corr and MAE table
    %     writecell({Vblnames{j}}, exportName, 'Sheet', j, 'Range', 'A1')
    %     writecell({'MAE'}, exportName, 'Sheet', j, 'Range', 'B2')
    %     writecell({'Model R^2 Values'}, exportName, 'Sheet', j, 'Range', 'C2');
    %     writecell({'Age log 10 P Values'}, exportName, 'Sheet', j, 'Range', 'D2')
    %     writecell({'Gender log 10 P Values'}, exportName, 'Sheet', j, 'Range', 'E2')
    %     writecell({'PC log 10 P Values'}, exportName, 'Sheet', j, 'Range', 'F2')
    %     writecell({'Age Only'}, exportName, 'Sheet', j, 'Range', 'A3');
    %     writecell({'Gender Only'}, exportName, 'Sheet', j, 'Range', 'A4');
    %     writecell({cornames{:}}', exportName, 'Sheet', j, 'Range', 'A5')
    %     writecell({'Samplenames'}, exportName, 'Sheet', j, 'Range', 'G2');
    %
end

%% export
% 
% writecell({'Average MAE LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 1)
% writecell({'Average R^2 Values LOOCV'}, exportName,'Range', 'A1', 'Sheet', 2);
% writecell({'Average Age Coefficient Estimates LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 3)
% writecell({'Average Gender Coefficient Estimates LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 4)
% writecell({'Average PC Coefficient Estimates LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 5)
% writecell({'Average Age log 10 P Values LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 6)
% writecell({'Average Gender log 10 P Values LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 7)
% writecell({'Average PC log 10 P Values LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 8)
writecell({'STD MAE LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 1)
writecell({'STD R^2 Values LOOCV'}, exportName,'Range', 'A1', 'Sheet', 2);
writecell({'STD Age Coefficient Estimates LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 3)
writecell({'STD Gender Coefficient Estimates LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 4)
writecell({'STD PC Coefficient Estimates LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 5)
writecell({'STD Age log 10 P Values LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 6)
writecell({'STD Gender log 10 P Values LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 7)
writecell({'STD PC log 10 P Values LOOCV'}, exportName, 'Range', 'A1', 'Sheet', 8)
for s = 1:8
    writecell(Vblnames, exportName, 'Range', 'B2', 'Sheet', s)
    writematrix(cornames, exportName, 'Range', 'A3', 'Sheet', s)
end
% writematrix(averageMAE, exportName, 'Range', 'B3', 'Sheet', 1)
writematrix(stdMAE, exportName, 'Range', 'B3', 'Sheet', 1)
% writematrix(averageRSQ, exportName, 'Range', 'B3', 'Sheet', 3)
writematrix(stdRSQ, exportName, 'Range', 'B3', 'Sheet', 2)
% writematrix(averageAgeP, exportName, 'Range', 'B3', 'Sheet', 5)
writematrix(stdAgeP, exportName, 'Range', 'B3', 'Sheet', 3)
% writematrix(averageGenderP, exportName, 'Range', 'B3', 'Sheet', 7)
writematrix(stdGenderP, exportName, 'Range', 'B3', 'Sheet', 4)
% writematrix(averagePCP, exportName, 'Range', 'B3', 'Sheet', 9)
writematrix(stdPCP, exportName, 'Range', 'B3', 'Sheet', 5)
% writematrix(averageAgeE, exportName, 'Range', 'B3', 'Sheet', 11)
writematrix(stdAgeE, exportName, 'Range', 'B3', 'Sheet', 6)
% writematrix(averageGenderE, exportName, 'Range', 'B3', 'Sheet', 13)
writematrix(stdGenderE, exportName, 'Range', 'B3', 'Sheet', 7)
% writematrix(averagePCE, exportName, 'Range', 'B3', 'Sheet', 15)
writematrix(stdPCE, exportName, 'Range', 'B3', 'Sheet', 8)      
       
% writecell({'PC LOOCV'}, exportName, 'Range', 'A1')
% writecell({'MAE LOOCV'}, exportName, 'Range', 'B1')
% writecell({'R^2 Values LOOCV'}, exportName,'Range', 'C1');
% writecell({'Age log 10 P Values'}, exportName, 'Range', 'D1')
% writecell({'Gender log 10 P Values'}, exportName, 'Range', 'E1')
% writecell({'PC log 10 P Values'}, exportName, 'Range', 'F1')
% writecell({'Traits'}, exportName, 'Range', 'G1');
% writematrix(reportData, exportName, 'Range', 'A2')
% writematrix(reportVbls, exportName, 'Range', 'G2')

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