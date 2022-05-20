pause
clear; clc; close;
%% Export names
exportName = 'PC_Trait_AssociationModel_TMFmethyl_FonlySites.xlsx';

%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
methylT = readtable('MethMatrixCG-HumanBuccalFitness192_10-100.csv', 'VariableNamingRule', 'preserve');
methylData = table2array(methylT(:, 2:end));
methylSites = methylT.Properties.VariableNames(2:end);
methylnames = methylT.sampleID;
convertNames = readtable('FitnessStudyIDkey_BarcodeIDs-to-SubjectIDs.csv', 'VariableNamingRule', 'preserve');
TM = readtable('TraitsAndDates_MarcoTomBuccal88_042722mjt', 'VariableNamingRule', 'preserve');

%% Import traits and import the sample names we have methylation values for

T = readtable('Traits_ProsperProteome103_120821mjt.csv', 'VariableNamingRule', 'preserve');
Vblnames = T.Properties.VariableNames(2:end);
traitdata = table2array(T(:, 2:end));
traitnames = T.sampleID;

%% check sites NOW TRY USING ONLY THIS METHYLATION DATA FOR THE PCs --> problem to think about is number of sites being consistent
proteomeT = readtable('MethMatrixCG-HumanBloodFitness96manual_20-100.csv', 'VariableNamingRule', 'preserve');
proteomeSites = proteomeT.Properties.VariableNames(2:end);


% for i = 1:numel(methylSites)
%     
%     for p = 1:numel(proteomeSites)
%         if strcmp(proteomeSites{p}, methylSites{i}) == 1 
%             sameSites(i) = p;
%             continue
%         end
%     end
% end
load('sameSites');
fitOnlySites = sameSites(find(sameSites ~= 0));
fitOnlySitesWithinAll = find(sameSites ~= 0);
%% find tom and marcos sampleIDs

for TMid = 1:height(TM.sampleID)
    for methID = 1:height(methylnames)
        if isequal(TM.sampleID(TMid), methylnames(methID))
            TMmethylIDs(TMid, 1) = methID;
        end
    end
end

m = 1; t = 1;
for h = 1:height(TM.subject)
   if isequal(TM.subject(h), {'marco'})
       marcoSampleID(m, 1) = TM.sampleID(h);
       m = m + 1;
   else
       tomSampleID(t, 1) = TM.sampleID(h);
       t = t + 1;
   end
end

for k=1:numel(methylnames)
    for j = 1:numel(tomSampleID)
        if isequal(tomSampleID(j), methylnames(k))
            tomID(j, 1) = k;
            continue
        end
    end
    for l = 1:numel(marcoSampleID)
        if isequal(marcoSampleID(l), methylnames(k))
            marcoID(l, 1) = k;
            continue
        end
    end       
end

marcoMethylData = methylData(marcoID, fitOnlySitesWithinAll);
tomMethylData = methylData(tomID, fitOnlySitesWithinAll);

%% Convert naming because the metadata and methylation data used different naming schemes
convertIdx = zeros(numel(methylnames),1);

for ID=1:numel(methylnames)
    for Cid =1:numel(convertNames.barcode)
        if isequal(convertNames.barcode(Cid), methylnames(ID))
            convertIdx(ID) = convertNames.subjectId(Cid);
            continue
        end
    end
end

samples2convert = find(convertIdx ~= 0);
samplenames(samples2convert) = convertIdx(samples2convert);
samplenames = samplenames(samples2convert)';
methylData = methylData(convertIdx(samples2convert), :);
%% Compare names from metadata and methylation value names to make sure they are the same
traitIdx = zeros(length(samplenames),1);

for kit=1:numel(samplenames)
    if length(find(traitnames == samplenames(kit))) == 1
        traitIdx(kit) = find(traitnames == samplenames(kit));
    end
end
%% PCA

mId = find(traitIdx ~= 0);
[methylCompsPCA, NEWscores] = pca(methylData(mId,fitOnlySitesWithinAll));
Y = traitdata(traitIdx(traitIdx ~= 0), :);

%% plots

cornames(:,1) = 1:width(methylCompsPCA);

MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
    ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];

[methylCompsVStraits, pval] = corr(Y, NEWscores);

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
    mae = zeros(width(NEWscores), 1);
    rSq = zeros(width(NEWscores), 1);
    agePval = zeros(width(NEWscores), 1);
    pcPval= zeros(width(NEWscores), 1);
    genderPval = zeros(width(NEWscores), 1);
    ageCoef = zeros(width(NEWscores), 1);
    pcCoef= zeros(width(NEWscores), 1);
    genderCoef = zeros(width(NEWscores), 1);
    
    for k = 1:width(NEWscores)
        X = [age, gender, NEWscores(:,k)];
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
        sig = find(pcPval(:,1) < -2);
              
    end
    
    if numel(sig) >= 1 && convertCharsToStrings(trait(:)) ~= "age" && convertCharsToStrings(trait(:)) ~= "sex"
        for p = 1:numel(sig)
            sigPC = sig(p);
            reportData(n,:) = [sigPC, mae(sigPC), rSq(sigPC), ageCoef(sigPC), agePval(sigPC), ...
                genderCoef(sigPC), genderPval(sigPC), pcCoef(sigPC), pcPval(sigPC), n];
            reportVbls(n,:) = convertCharsToStrings({trait(:)});
            n = n + 1;
        end
        
        
    end
    
end

%% Export PC associated trait data

reportData = sortrows(reportData, 1);
reportVbls = reportVbls(reportData(:,10));

writecell({'PC'}, exportName, 'Range', 'A1')
writecell({'MAE'}, exportName, 'Range', 'B1')
writecell({'R^2 Values'}, exportName,'Range', 'C1');
writecell({'Age Coef Values'}, exportName, 'Range', 'D1')
writecell({'Age log 10 P Values'}, exportName, 'Range', 'E1')
writecell({'Gender Coef Values'}, exportName, 'Range', 'F1')
writecell({'Gender log 10 P Values'}, exportName, 'Range', 'G1')
writecell({'PC Coef Values'}, exportName, 'Range', 'H1')
writecell({'PC log 10 P Values'}, exportName, 'Range', 'I1')

writecell({'Traits'}, exportName, 'Range', 'J1');
writematrix(reportData(:,1:9), exportName, 'Range', 'A2')
writematrix(reportVbls, exportName, 'Range', 'J2');
writecell({'=IF(A2=A1,K1,K1+1)'}, exportName, 'Range', 'K2');

%% Predict Tom and Marcos Scores

tomScores = tomMethylData * methylCompPCA;
marcoScores = marcoMethylData * methylCompPCA;

sigPCs = 2; % This needs to be determined

for p = 1:numel(sigPCs)
    fig = figure(1);
    x = [NEWscores(:,sigPCs(p)); tomScores(:,sigPCs(p)); marcoScores(:,sigPCs(p))];
    l1 = repmat({'Fitness Study'}, height(NEWscores), 1);
    l2 = repmat({'Tom'}, height(tomScores), 1);
    l3 = repmat({'Marco'}, height(marcoScores), 1);
    l = [l1; l2; l3];
    boxplot(x, l);
    
    xl = sprintf("PC %.0f", sigPCs(p));
    xlabel(xl);
    ylabel("Component Value");
    savename = sprintf('boxplots_PC%d.png', sigPCs(p));
    
    saveas(fig, savename);
    close all;
    clear fig;
end