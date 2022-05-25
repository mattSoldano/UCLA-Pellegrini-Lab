pause
clear; clc; close;
%% Export names
exportName = 'SitePredict_Buccal.xlsx';

%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
buccalT = readtable('MethMatrixCG-HumanBuccalFitness192_10-100.csv', 'VariableNamingRule', 'preserve');
buccalData = table2array(buccalT(:, 2:end));
buccalSites = buccalT.Properties.VariableNames(2:end);
samplenames = buccalT.sampleID;

%% Import traits and import the sample names we have methylation values for

T = readtable('Traits_HumanBuccalFitness104x173_040422mjt.csv', 'VariableNamingRule', 'preserve');
TM = readtable('TraitsAndDates_MarcoTomBuccal88_042722mjt', 'VariableNamingRule', 'preserve');

Vblnames = T.Properties.VariableNames(2:end);
traitdata = table2array(T(:, 2:end));
traitnames = T.subjectID;
%% Compare names from metadata and methylation value names to make sure they are the same
traitIdx = zeros(length(samplenames),1);

for kit=1:numel(traitnames)
    if length(find(strcmp(traitnames(kit), samplenames))) == 1
        traitIdx(kit) = find(strcmp(traitnames(kit), samplenames));
    end
end
%% Sort data

fitIDs = find(traitIdx ~= 0);
% fitIDs = fitness trait sample indices
mS = traitIdx(fitIDs);
% mS --> IE the methylation sample indices in order of
% the trait data IDs
% --> confirmed that samplenames(mS) = traitnames(fitIDs)
% --> IE the methylation sample indices and trait sample indices are aligned

sortedBuccalData = buccalData(mS,:);
Y = traitdata(fitIDs, :);

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
    
    for k = 1:width(sortedBuccalData)
        X = [age(l3d.(trait)), gender(l3d.(trait)), sortedBuccalData(l3d.(trait),k)];
        mdl = fitlm(X , Y(l3d.(trait),j));
        sitePval = log10(mdl.Coefficients.pValue(4));
        
        if isinf(sitePval)
            disp('hi');
        end
        
        if sitePval < -5 && sitePval > -inf
            % level of significance = 10^-5
            site = buccalSites{k};
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

cd BuccalSiteSpecificresults
save('BuccalSiteSpecific_Trait_Associations');

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

marcoMethylData = buccalData(marcoID, :);
tomMethylData = buccalData(tomID, :);
%% Predict Tom and Marcos Scores

sortData = sortrows(reportData, 9);
[u, ia, ic] = unique(sortData(:,1));
uniqueSigSites = u((ic),1);

for nu = 1:numel(uniqueSigSites)
    if numel(find(strcmp(uniqueSigSites(nu), buccalSites))) == 1
        idxUniqueSigSites(nu) = find(strcmp(uniqueSigSites(nu), buccalSites));
    end
end

for p = 1:10 % top 10 for now
        
    uniqueMarcoData = unique(marcoMethylData(:,idxUniqueSigSites(p)));
    uniqueTomData = unique(tomMethylData(:,idxUniqueSigSites(p)));
    uniqueBuccalData = unique(sortedBuccalData(:,idxUniqueSigSites(p)));
    
    x = [uniqueMarcoData; uniqueTomData; uniqueBuccalData];
    
    l1 = repmat({'1 '}, height(uniqueMarcoData), 1);
    % Marco
    l2 = repmat({'2 '}, height(uniqueTomData), 1);
    % tom    
    l3 = repmat({'Fitness Study'}, height(uniqueBuccalData), 1);

    l = [l1; l2; l3];
    
    siteplot = uniqueSigSites{p};
    
    cd /Users/Sully/Pellegrini/FitnessData
    fig = figure(1);
    if numel(find(x == 0)) > 1
        boxplot(x, l);
        savename = sprintf('Buccal_Boxplot_%s.png', siteplot);
    else
        violinplot(x, l);
        savename = sprintf('Buccal_Violinplot_%s.png', siteplot);
    end

    siteplot = regexprep(siteplot, '_', ':');
    xlabel(siteplot, 'FontSize', 15);
    ylabel("Methylation Value", 'FontSize', 15);
    
    cd BuccalSiteSpecificresults
    saveas(fig, savename);
    close all;
    clear fig;
end