pause
clear; clc; close;
%% Export names
exportName = 'tripleCheck_SitePredict.xlsx';

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
    
    for k = 1:width(fitMethylData)
        X = [age, gender, fitMethylData(:,k)];
        mdl = fitlm(X , Y(:,j));
        sitePval = log10(mdl.Coefficients.pValue(4));
        
        if sitePval < -5
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
allData=[reportData(:,1:10), {PresentreportVbls{:}}'];
save('SiteSpecificAssocationsAll');

writecell({'Site'}, exportName, 'Range', 'A1')
writecell({'MAE'}, exportName, 'Range', 'B1')
writecell({'R^2 Values'}, exportName,'Range', 'C1');
writecell({'Age Coef Values'}, exportName, 'Range', 'D1')
writecell({'Age log 10 P Values'}, exportName, 'Range', 'E1')
writecell({'Gender Coef Values'}, exportName, 'Range', 'F1')
writecell({'Gender log 10 P Values'}, exportName, 'Range', 'G1')
writecell({'PC Coef Values'}, exportName, 'Range', 'H1')
writecell({'PC log 10 P Values'}, exportName, 'Range', 'I1')
writecell({'Intercept'}, exportName, 'Range', 'J1')
writecell({'Traits'}, exportName, 'Range', 'K1');
writecell(allData, exportName, 'Range', 'A2');