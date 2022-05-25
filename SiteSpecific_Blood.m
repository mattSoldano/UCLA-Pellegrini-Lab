pause
clear; clc; close;
%% Export names
exportName = 'SitePredict_Blood.xlsx';

%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
methylT = readtable('MethMatrixCG-HumanBuccalFitness192_10-100.csv', 'VariableNamingRule', 'preserve');
methylData = table2array(methylT(:, 2:end));
methylSites = methylT.Properties.VariableNames(2:end);
samplenames = methylT.sampleID;

%% Import traits and import the sample names we have methylation values for

T = readtable('Traits_HumanBuccalFitness104x173_040422mjt.csv', 'VariableNamingRule', 'preserve');
TM = readtable('TraitsAndDates_MarcoTomBuccal88_042722mjt', 'VariableNamingRule', 'preserve');
convertNames = readtable('FitnessStudyIDkey_BarcodeIDs-to-SubjectIDs.csv', 'VariableNamingRule', 'preserve');

Vblnames = T.Properties.VariableNames(2:end);
traitdata = table2array(T(:, 2:end));
traitnames = T.subjectID;

%% Read in blood data
bloodT = readtable('MethMatrixCG-HumanBloodFitness96manual_20-100.csv', 'VariableNamingRule', 'preserve');
bloodSites = bloodT.Properties.VariableNames(2:end);
bloodData = table2array(bloodT(:, 2:end));
bloodnames = bloodT.sampleId;

load('sameSites');
bloodSameSitesIdx = sameSites(find(sameSites ~= 0));
% for blood data
buccalSameSitesIdx = find(sameSites ~= 0);
% for buccal data

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

SortedBloodData = bloodData(bloodIdx, bloodSameSitesIdx);
Y = Y(traitIdx(bloodIdx), :);
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
    
    l3d.(trait) = find(~isoutlier(Y(:,j)));
    % less than 3 standard deviations indicies (remove outliers)
    
    for k = 1:width(SortedBloodData)
        X = [age(l3d.(trait)), gender(l3d.(trait)), SortedBloodData(l3d.(trait),k)];
        mdl = fitlm(X , Y(l3d.(trait),j));
        sitePval = log10(mdl.Coefficients.pValue(4));
        
        if sitePval < -3 && sitePval > -inf
            % level of significance = 10^-5
            site = methylSites{k};
            mae = sqrt(mdl.MSE);
            rSq = mdl.Rsquared.Ordinary;
            agePval = log10(mdl.Coefficients.pValue(2));
            genderPval = log10(mdl.Coefficients.pValue(3));
            intercept = mdl.Coefficients.Estimate(1);
            ageCoef = mdl.Coefficients.Estimate(2);
            genderCoef = mdl.Coefficients.Estimate(3);
            siteCoef = mdl.Coefficients.Estimate(4);
            
            reportData(n,:) = [{site}, mae, rSq, ageCoef, agePval, ...
                genderCoef, genderPval, siteCoef, sitePval, ...
                intercept,n];
            reportVbls(n,:) = convertCharsToStrings({trait(:)});
            
            beta(n,:) = [ageCoef,genderCoef,siteCoef];
            predictY(n,:) = [{trait}; {site}; (X * beta(n,:)') + intercept];
           n = n + 1;
        end
        
    end
    
end

%% Export PC associated trait data

PresentreportVbls = reportVbls(cell2mat(reportData(:,11)));
allData=[num2cell(reportData(:,1:10)), {PresentreportVbls{:}}'];

cd BloodSiteSpecificresults
save('BloodSiteSpecific_Trait_Associations');

writecell({'PC'}, exportName, 'Range', 'A1')
writecell({'MAE'}, exportName, 'Range', 'B1')
writecell({'Model R^2'}, exportName,'Range', 'C1');
writecell({'Age Coef'}, exportName, 'Range', 'D1')
writecell({'Age log 10 P Values'}, exportName, 'Range', 'E1')
writecell({'Gender Coef'}, exportName, 'Range', 'F1')
writecell({'Gender log 10 P Values'}, exportName, 'Range', 'G1')
writecell({'Site Coef'}, exportName, 'Range', 'H1')
writecell({'Site log 10 P Values'}, exportName, 'Range', 'I1')
writecell({'Intercept'}, exportName, 'Range', 'J1')
writecell({'Traits'}, exportName, 'Range', 'K1');
writecell(reportData(:,1:10), exportName, 'Range', 'A2')
writematrix(reportVbls, exportName, 'Range', 'K2');
writecell({'=IF(A2=A1,L1,L1+1)'}, exportName, 'Range', 'L2');

%% find tom and marcos sampleIDs

for TMid = 1:height(TM.sampleID)
    for methID = 1:height(samplenames)
        if isequal(TM.sampleID(TMid), samplenames(methID))
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

for k=1:numel(samplenames)
    for j = 1:numel(tomSampleID)
        if isequal(tomSampleID(j), samplenames(k))
            tomID(j, 1) = k;
            continue
        end
    end
    for l = 1:numel(marcoSampleID)
        if isequal(marcoSampleID(l), samplenames(k))
            marcoID(l, 1) = k;
            continue
        end
    end       
end

marcoMethylData = methylData(marcoID, buccalSameSitesIdx);
tomMethylData = methylData(tomID, buccalSameSitesIdx);

%% Predict Tom and Marcos Scores

% uniqueSigSites = unique(reportData(:,1));
% 
% for nu = 1:numel(uniqueSigSites)
%     if numel(find(strcmp(uniqueSigSites(nu), methylSites))) == 1
%         idxUniqueSigSites(nu) = find(strcmp(uniqueSigSites(nu), methylSites));
%     end
% end

% Temporary for seeing violin plots of sites other than neutrophil sites

sortData = sortrows(reportData, 9);
top10id = [1,2,3,80,81,136,166,243,339,361];
tempUniqSigSites = sortData(top10id,1);

for nu = 1:numel(tempUniqSigSites)
    if numel(find(strcmp(tempUniqSigSites(nu), methylSites))) == 1
        tempIdxUniqSigSites(nu) = find(strcmp(tempUniqSigSites(nu), methylSites));
    end
end

for p = 1:numel(tempUniqSigSites) % Top 10 interesting traits (not all neutrophils)
    fig = figure(1);
    
    uniqueMarcoData = unique(marcoMethylData(:,top10id(p)));
    uniqueTomData = unique(tomMethylData(:,top10id(p)));
    uniqueBloodData = unique(SortedBloodData(:,top10id(p)));
    
    x = [uniqueMarcoData; uniqueTomData; uniqueBloodData];
    
    l1 = repmat({'1 '}, height(uniqueMarcoData), 1);
    % Marco
    l2 = repmat({'2 '}, height(uniqueTomData), 1);
    % tom    
    l3 = repmat({'Fitness Study'}, height(uniqueBloodData), 1);

    l = [l1; l2; l3];
    
    siteplot = tempUniqSigSites{p};
    
    cd /Users/Sully/Pellegrini/FitnessData        
    if numel(find(x == 0)) > 1
        boxplot(x, l);
        savename = sprintf('Blood_Boxplot_%s.png', siteplot);
    else
        violinplot(x, l);
        savename = sprintf('Blood_Violinplot_%s.png', siteplot);
    end
           
    siteplot = regexprep(siteplot, '_', ':');
    xlabel(siteplot, 'FontSize', 15);
    ylabel("Methylation Value", 'FontSize', 15);
    
    cd BloodSiteSpecificresults
    saveas(fig, savename);
    close all;
    clear fig;
end