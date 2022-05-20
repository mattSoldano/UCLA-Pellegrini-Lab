pause
clear; clc; close;
%% Export names
exportName = 'PC_Trait_AssociationModel_allMethylData.xlsx';

%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
methylT = readtable('MethMatrixCG-HumanBuccalFitness192_10-100.csv', 'VariableNamingRule', 'preserve');
methylData = table2array(methylT(:, 2:end));
methylSites = methylT.Properties.VariableNames(2:end);
samplenames = methylT.sampleID;

%% Import traits and import the sample names we have methylation values for

T = readtable('Traits_HumanBuccalFitness104x173_040422mjt.csv', 'VariableNamingRule', 'preserve');
TM = readtable('TraitsAndDates_MarcoTomBuccal88_042722mjt', 'VariableNamingRule', 'preserve');

Vblnames = T.Properties.VariableNames(2:end);
traitdata = table2array(T(:, 2:end));
traitnames = T.subjectID;
%% find tom and marcos sampleIDs

m = 1; t = 1;
for h = 1:height(TM.subject)
   if strcmp(TM.subject(h), {'marco'})
       marcoSampleID(m, 1) = TM.sampleID(h);
       m = m + 1;
   elseif strcmp(TM.subject(h), {'tom'})
       tomSampleID(t, 1) = TM.sampleID(h);
       t = t + 1;
   end
end

for k=1:numel(samplenames)
    for j = 1:numel(tomSampleID)
        if strcmp(tomSampleID(j), samplenames(k))
            tomID(j, 1) = k;
            continue
        end
    end
    for l = 1:numel(marcoSampleID)
        if strcmp(marcoSampleID(l), samplenames(k))
            marcoID(l, 1) = k;
            continue
        end
    end       
end

marcoMethylData = methylData(marcoID, :);
tomMethylData = methylData(tomID, :);
%% Compare names from metadata and methylation value names to make sure they are the same
traitIdx = zeros(length(samplenames),1);

for kit=1:numel(traitnames)
    if length(find(strcmp(traitnames(kit), samplenames))) == 1
        traitIdx(kit) = find(strcmp(traitnames(kit), samplenames));
    end
end
%% PCA

fitIDs = find(traitIdx ~= 0);
% fitIDs = fitness trait sample indices
mS = traitIdx(fitIDs);
% mS --> IE the methylation sample indices in order of
% the trait data IDs
% --> confirmed that samplenames(mS) = traitnames(fitIDs)
    % --> IE the methylation sample indices and trait sample indices are aligned

fitMethylData = methylData(mS,:);
Y = traitdata(fitIDs, :);
marcoTomFitMethylData = [fitMethylData; tomMethylData; marcoMethylData];
[methylLoadings, PCscores, latent, tsquared, explained, mu] = pca(marcoTomFitMethylData);
PCscoresRecover = (marcoTomFitMethylData - mu) * methylLoadings;

%% plots

cornames(:,1) = 1:width(methylLoadings);
%% Data

age = Y(:,3);
gender = Y(:,2);
n = 1;
for j = 1:width(Vblnames)
    j
    trait = Vblnames{j};
    trait = regexprep(trait, ' ', '_');
    trait = regexprep(trait, '-', '_');
    trait = regexprep(trait, '%', 'percent');
    trait = regexprep(trait, '#', 'number');
    trait = regexprep(trait, '\.', '_');
    mae = zeros(width(PCscores), 1);
    rSq = zeros(width(PCscores), 1);
    agePval = zeros(width(PCscores), 1);
    pcPval= zeros(width(PCscores), 1);
    genderPval = zeros(width(PCscores), 1);
    ageCoef = zeros(width(PCscores), 1);
    pcCoef= zeros(width(PCscores), 1);
    genderCoef = zeros(width(PCscores), 1);
    intercept = zeros(width(methylData), 1);
    
    for k = 1:width(PCscores)
        X = [age, gender, PCscores(1:numel(mS),k)];
        mdl = fitlm(X , Y(:,j));
        mdlCoefficients.(trait) = mdl.Coefficients;
        mae(k) = sqrt(mdl.MSE);
        rSq(k) = mdl.Rsquared.Ordinary;
        agePval(k) = log10(mdlCoefficients.(trait).pValue(2));
        genderPval(k) = log10(mdlCoefficients.(trait).pValue(3));
        pcPval(k) = log10(mdlCoefficients.(trait).pValue(4));
        ageCoef(k) = mdlCoefficients.(trait).Estimate(2);
        genderCoef(k) = mdlCoefficients.(trait).Estimate(3);
        pcCoef(k) = mdlCoefficients.(trait).Estimate(4);
        intercept(k) = mdlCoefficients.(trait).Estimate(1);
        sig = find(pcPval(:,1) < -5);
        % level of significance ^
              
    end
    
    if numel(sig) >= 1 && convertCharsToStrings(trait(:)) ~= "age" && convertCharsToStrings(trait(:)) ~= "sex"
        for p = 1:numel(sig)
            sigPC(n) = sig(p);
            reportData(currentn,:) = [sigSite(currentn), mae(sigSite(currentn)), ...
                rSq(sigSite(currentn)), ageCoef(sigSite(currentn)), agePval(sigSite(currentn)), ...
                genderCoef(sigSite(currentn)), genderPval(sigSite(currentn)), ...
                siteCoef(sigSite(currentn)), sitePval(sigSite(currentn)), ...
                intercept(sigSite(currentn)), currentn];
            reportVbls(n,:) = convertCharsToStrings({trait(:)});
            n = n + 1;            
        end
    end
end

%% Export PC associated trait data

reportData = sortrows(reportData, 1);
reportVbls = reportVbls(reportData(:,11));

z = 1;
shade = zeros(height(reportData),1);
shade(1) = 1;
for h = 1:height(reportData)-1
    if reportData(h+1,1) > reportData(h,1)
        shade(h+1) =  z + 1;
        z = z + 1;
    else
        shade(h+1) = z;
    end
end

writecell({'PC'}, exportName, 'Range', 'A1')
writecell({'MAE'}, exportName, 'Range', 'B1')
writecell({'R^2 Values'}, exportName,'Range', 'C1');
writecell({'Age Coef Values'}, exportName, 'Range', 'D1')
writecell({'Age log 10 P Values'}, exportName, 'Range', 'E1')
writecell({'Gender Coef Values'}, exportName, 'Range', 'F1')
writecell({'Gender log 10 P Values'}, exportName, 'Range', 'G1')
writecell({'PC Coef Values'}, exportName, 'Range', 'H1')
writecell({'PC log 10 P Values'}, exportName, 'Range', 'I1')
writecell({'Intercept'}, exportName,'Range', 'J1');

writecell({'Traits'}, exportName, 'Range', 'K1');
writematrix(reportData(:,1:10), exportName, 'Range', 'A2')
writematrix(reportVbls, exportName, 'Range', 'K2');
writematrix(shade, exportName, 'Range', 'L2');

%% Predict Tom and Marcos Scores

% tomScores = tomMethylData * methylLoadings;
% marcoScores = marcoMethylData * methylLoadings;
% 
% tomScoresIDs = numel(mS):numel(mS) + numel(tomID);
% marcoScoresIDs = numel(mS) + (numel(tomID)): numel(mS) + numel(tomID)+ numel(marcoID);
% 
% uniqueSigPCs = unique(sigPC);
% uniqueSigPCs = [9, 19, 33, 44, 80, 104];
% 
% for p = 1:numel(uniqueSigPCs)
%     fig = figure('visible','off');
% %     x = [PCscores(:,uniqueSigPCs(p)); tomScores(:,uniqueSigPCs(p)); marcoScores(:,uniqueSigPCs(p))];
%     x = [PCscores(1:numel(mS),uniqueSigPCs(p)); PCscores(tomScoresIDs,uniqueSigPCs(p)); ...
%         PCscores(marcoScoresIDs,uniqueSigPCs(p))];
%     l1 = repmat({'Fitness Study'}, numel(mS), 1);
%     l2 = repmat({'1 '}, numel(tomScoresIDs), 1);
%     l3 = repmat({'2 '}, numel(marcoScoresIDs), 1);
%     l = [l1; l2; l3];
%     violinplot(x, l);
%     
%     xl = sprintf("PC %.0f", uniqueSigPCs(p));
%     xlabel(xl);
%     ylabel("Component Value");
%     savename = sprintf('a_Violinplots_AllDataPCs_PC%d.png', uniqueSigPCs(p));
%     
%     saveas(fig, savename);
%     close all;
%     clear fig;
% end
% 
% %     fig = figure('visible','off');
