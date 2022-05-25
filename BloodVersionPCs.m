pause
clear; clc; close;
%% Export names
exportName = 'PC_Trait_AssociationModel_Blood.xlsx';

%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
buccalT = readtable('MethMatrixCG-HumanBuccalFitness192_10-100.csv', 'VariableNamingRule', 'preserve');
buccalData = table2array(buccalT(:, 2:end));
buccalSites = buccalT.Properties.VariableNames(2:end);
buccalnames = buccalT.sampleID;
convertNames = readtable('FitnessStudyIDkey_BarcodeIDs-to-SubjectIDs.csv', 'VariableNamingRule', 'preserve');
TM = readtable('TraitsAndDates_MarcoTomBuccal88_042722mjt', 'VariableNamingRule', 'preserve');

%% Import traits and import the sample names we have methylation values for

T = readtable('Traits_HumanBuccalFitness104x173_040422mjt.csv', 'VariableNamingRule', 'preserve');
Vblnames = T.Properties.VariableNames(2:end);
traitdata = table2array(T(:, 2:end));
traitnames = T.subjectID;

%% check sites NOW TRY USING ONLY THIS METHYLATION DATA FOR THE PCs --> problem to think about is number of sites being consistent
bloodT = readtable('MethMatrixCG-HumanBloodFitness96manual_20-100.csv', 'VariableNamingRule', 'preserve');
bloodSites = bloodT.Properties.VariableNames(2:end);
bloodData = table2array(bloodT(:, 2:end));
bloodnames = bloodT.sampleId;

% for i = 1:numel(buccalSites)
%     
%     for p = 1:numel(bloodSites)
%         if strcmp(bloodSites{p}, buccalSites{i}) == 1 
%             sameSites(i) = p;
%             continue
%         end
%     end
% end
load('sameSites');
bloodSameSitesIdx = sameSites(find(sameSites ~= 0));
% for blood data
buccalSameSitesIdx = find(sameSites ~= 0);
% for buccal data
%% find tom and marcos sampleIDs

for TMid = 1:height(TM.sampleID)
    for methID = 1:height(buccalnames)
        if isequal(TM.sampleID(TMid), buccalnames(methID))
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

for k=1:numel(buccalnames)
    for j = 1:numel(tomSampleID)
        if isequal(tomSampleID(j), buccalnames(k))
            tomID(j, 1) = k;
            continue
        end
    end
    for l = 1:numel(marcoSampleID)
        if isequal(marcoSampleID(l), buccalnames(k))
            marcoID(l, 1) = k;
            continue
        end
    end       
end

marcoMethylData = buccalData(marcoID, buccalSameSitesIdx);
tomMethylData = buccalData(tomID, buccalSameSitesIdx);

%% Convert naming because the metadata and methylation data used different naming schemes
convert = zeros(numel(traitnames),1);

for ID=1:numel(traitnames)
    for Cid =1:numel(convertNames.barcode)
        if isequal(convertNames.barcode(Cid), traitnames(ID))
            convert(ID) = convertNames.subjectId(Cid);
            continue
        end
    end
end

samples2convert = find(convert ~= 0);
traitnames = convert(samples2convert);
% convert = id's at converted indicies
% samples2convert = indicies converted
% --> traitnames and Y are aligned
Y = traitdata(samples2convert, :);
%% Compare names from metadata and methylation value names to make sure they are the same
traitIdx = zeros(length(bloodnames),1);

for kit=1:numel(bloodnames)
    if length(find(traitnames == bloodnames(kit))) == 1
        traitIdx(kit) = find(traitnames == bloodnames(kit));
    end
end

bloodIdx = find(traitIdx ~= 0);
newBloodnames = bloodnames(bloodIdx);
newTraitnames = traitnames(traitIdx(bloodIdx));

% traits(bloodIdx(traitIdx)) correspond to blooddata(traitIdx)

bloodData = bloodData(bloodIdx, bloodSameSitesIdx);
Y = Y(traitIdx(bloodIdx), :);
%% PCA

[methylCompsPCA, BloodScores] = pca(bloodData);

%% Data

age = Y(:,3);
gender = Y(:,2);
n = 1;
for j = 1:width(Vblnames)
    trait = Vblnames{j};
    fprintf('%.0f:%s\n', j, trait);
    trait = regexprep(trait, ' ', '_');
    trait = regexprep(trait, '-', '_');
    trait = regexprep(trait, '%', 'percent');
    trait = regexprep(trait, '#', 'number');
    trait = regexprep(trait, '\.', '_');
    
    if convertCharsToStrings(trait(:)) == "age" ...
            || convertCharsToStrings(trait(:)) == "sex"
        continue
    end
    
    if convertCharsToStrings(trait(:)) == "group"
        l3d.(trait) = 1:numel(Y(:,j));
    else 
        l3d.(trait) = find(~isoutlier(Y(:,j)));
    end
    % less than 3 standard deviations indicies (remove outliers)
    
    for k = 1:width(BloodScores)
        X = [age(l3d.(trait)), gender(l3d.(trait)), BloodScores(l3d.(trait),k)];
        mdl = fitlm(X , Y(l3d.(trait),j));
        PC_Pval = log10(mdl.Coefficients.pValue(4));
        
        if PC_Pval < -3 && PC_Pval > -inf
            % level of significance = 10^-5
            PC = k;
            mae = sqrt(mdl.MSE);
            rSq = mdl.Rsquared.Ordinary;
            agePval = log10(mdl.Coefficients.pValue(2));
            genderPval = log10(mdl.Coefficients.pValue(3));
            intercept = mdl.Coefficients.Estimate(1);
            ageCoef = mdl.Coefficients.Estimate(2);
            genderCoef = mdl.Coefficients.Estimate(3);
            PC_Coef = mdl.Coefficients.Estimate(4);
            
            reportData(n,:) = [PC, mae, rSq, ageCoef, agePval, ...
                genderCoef, genderPval, PC_Coef, PC_Pval, ...
                intercept,n];
            reportVbls(n,:) = convertCharsToStrings({trait(:)});
            
            beta(n,:) = [ageCoef,genderCoef,PC_Coef];
            predictY(n,:) = [{trait}; PC; (X * beta(n,:)') + intercept];
            n = n + 1;
        end
        
    end
    
end

%% Export PC associated trait data
PresentreportVbls = reportVbls(reportData(:,11));
allData=[num2cell(reportData(:,1:10)), {PresentreportVbls{:}}'];

cd BloodPCresults
save('bloodMethylPC_Trait_Associations');

writecell({'PC'}, exportName, 'Range', 'A1')
writecell({'MAE'}, exportName, 'Range', 'B1')
writecell({'Model R^2'}, exportName,'Range', 'C1');
writecell({'Age Coef'}, exportName, 'Range', 'D1')
writecell({'Age log 10 P Values'}, exportName, 'Range', 'E1')
writecell({'Gender Coef'}, exportName, 'Range', 'F1')
writecell({'Gender log 10 P Values'}, exportName, 'Range', 'G1')
writecell({'PC Coef'}, exportName, 'Range', 'H1')
writecell({'PC log 10 P Values'}, exportName, 'Range', 'I1')
writecell({'Intercept'}, exportName, 'Range', 'J1')
writecell({'Traits'}, exportName, 'Range', 'K1');
writematrix(reportData(:,1:10), exportName, 'Range', 'A2')
writematrix(reportVbls, exportName, 'Range', 'K2');
writecell({'=IF(A2=A1,L1,L1+1)'}, exportName, 'Range', 'L2');

%% Predict Tom and Marcos Scores

tomScores = tomMethylData * methylCompsPCA;
marcoScores = marcoMethylData * methylCompsPCA;

uniqueSigPCs = unique(reportData(:,1));

for p = 1:numel(uniqueSigPCs)
    fig = figure(1);
    x = [marcoScores(:,uniqueSigPCs(p)); tomScores(:,uniqueSigPCs(p)); BloodScores(:,uniqueSigPCs(p))];
    
    l1 = repmat({'1 '}, height(marcoScores), 1);
    % Marco    
    l2 = repmat({'2 '}, height(tomScores), 1);
    % tom
    l3 = repmat({'Fitness Study'}, height(BloodScores), 1);
    
    l = [l1; l2; l3];
    cd ..
    violinplot(x, l);
    
    xl = sprintf("PC %.0f", uniqueSigPCs(p));
    xlabel(xl);
    ylabel("Component Value");
    savename = sprintf('Blood_Violinplots_PC%d.png', uniqueSigPCs(p));
    
    cd BloodPCresults
    saveas(fig, savename);
    close all;
    clear fig;
end