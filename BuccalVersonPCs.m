pause
clear; clc; close;
%% Export names
exportName = 'PC_Trait_AssociationModel_Buccal.xlsx';

%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
buccalT = readtable('MethMatrixCG-HumanBuccalFitness192_10-100.csv', 'VariableNamingRule', 'preserve');
buccalData = table2array(buccalT(:, 2:end));
buccalSites = buccalT.Properties.VariableNames(2:end);
buccalnames = buccalT.sampleID;
TM = readtable('TraitsAndDates_MarcoTomBuccal88_042722mjt', 'VariableNamingRule', 'preserve');

%% Import traits and import the sample names we have methylation values for

T = readtable('Traits_HumanBuccalFitness104x173_040422mjt.csv', 'VariableNamingRule', 'preserve');
Vblnames = T.Properties.VariableNames(2:end);
traitdata = table2array(T(:, 2:end));
traitnames = T.subjectID;
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

marcoMethylData = buccalData(marcoID, :);
tomMethylData = buccalData(tomID, :);

%% Compare names from metadata and methylation value names to make sure they are the same
traitIdx = zeros(length(buccalnames),1);

for kit=1:numel(buccalnames)
    for buccali=1:numel(traitnames)
        if length(find(isequal(traitnames{buccali}, buccalnames{kit}))) == 1
            traitIdx(kit) = buccali;
        end
    end
end

buccalidx = find(traitIdx ~= 0);

buccalnames = buccalnames(buccalidx);
buccalData = buccalData(buccalidx, :);

traitnames = traitnames(traitIdx(buccalidx));
Y = traitdata(traitIdx(buccalidx), :);

if isequal(traitnames,buccalnames) == false
    warning('sample and trait IDs don''t match');
end
% traits(traitIdx(buccalidx)) correspond to buccaldata(buccalidx)

%% PCA

[methylCompsPCA, BuccalScores] = pca(buccalData);
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
    
    for k = 1:width(BuccalScores)
        X = [age(l3d.(trait)), gender(l3d.(trait)), BuccalScores(l3d.(trait),k)];
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

cd BuccalPCresults
save('buccalMethylPC_Trait_Associations');

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
    x = [marcoScores(:,uniqueSigPCs(p)); tomScores(:,uniqueSigPCs(p)); BuccalScores(:,uniqueSigPCs(p))];
    
    l1 = repmat({'1 '}, height(marcoScores), 1);
    % Marco    
    l2 = repmat({'2 '}, height(tomScores), 1);
    % tom
    l3 = repmat({'Fitness Study'}, height(BuccalScores), 1);
    
    l = [l1; l2; l3];
    cd ..
    violinplot(x, l);
    
    xl = sprintf("PC %.0f", uniqueSigPCs(p));
    xlabel(xl);
    ylabel("Component Value");
    savename = sprintf('buccal_Violinplots_PC%d.png', uniqueSigPCs(p));
    
    cd BuccalPCresults
    saveas(fig, savename);
    close all;
    clear fig;
end