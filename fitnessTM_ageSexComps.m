pause
clear; clc; close;
%% Export names
exportName = 'Fitness_TomMarco_PC_Age_Gender_predictTraits3.xlsx';

%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
methylT = readtable('MethMatrixCG-HumanBuccalFitness192_10-100.csv', 'VariableNamingRule', 'preserve');
methylData = table2array(methylT(:, 2:end));
methylSites = methylT.Properties.VariableNames(2:end);
samplenames = methylT.sampleID;
%% Import traits and import the sample names we have methylation values for

T = readtable('Traits_HumanBuccalFitness104x173_040422mjt.csv', 'VariableNamingRule', 'preserve');
TM = readtable('TraitsAndDates_MarcoTomBuccal88_042722mjt', 'VariableNamingRule', 'preserve');
TMage = TM.age; TMsex = TM.sex;
g3 = find(TMage == 33);
g4 = find(TMage ~= 33);
% TMgroup(g3, 1) = 3; TMgroup(g4, 1) = 4;
% Marco is group 3, Tom is group 4
Vblnames = T.Properties.VariableNames(2:end);
traitdata = table2array(T(:, 2:end));
TMdata = NaN(height(TMage), width(traitdata));
% TMdata(:,[2:3]) = [TMsex, TMage];
traitdata = [traitdata; TMdata];
traitnames = {T.subjectID{:}}';%, TM.sampleID{:}}';

%% Convert names

%% Compare names from metadata and methylation value names to make sure they are the same
traitIdx = zeros(length(samplenames),1);

for kit=1:numel(samplenames)
    if length(find(strcmp(traitnames,samplenames{kit}))) == 1
        traitIdx(kit) = find(strcmp(traitnames,samplenames{kit}));
    end
end
%% PCA

mId = find(traitIdx ~= 0);
[methylCompsPCA, scores] = pca(methylData(mId,:));
Y = traitdata(traitIdx ~= 0, :);


%% plots

cornames(:,1) = 1:width(methylCompsPCA);

MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
    ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];

[methylCompsVStraits, pval] = corr(Y, scores);

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
    mae = zeros(width(scores), 1);
    rSq = zeros(width(scores), 1);
    agePval = zeros(width(scores), 1);
    pcPval= zeros(width(scores), 1);
    genderPval = zeros(width(scores), 1);
    ageCoef = zeros(width(scores), 1);
    pcCoef= zeros(width(scores), 1);
    genderCoef = zeros(width(scores), 1);
    
    for k = 1:width(scores)
        X = [age, gender, scores(:,k)];
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
            reportData(n,:) = [sigPC, mae(sigPC), rSq(sigPC), agePval(sigPC), ageCoef(sigPC), ...
                genderPval(sigPC), genderCoef(sigPC), pcPval(sigPC), pcCoef(sigPC), n];
            reportVbls(n,:) = convertCharsToStrings({trait(:)});
            n = n + 1;
        end
        
        
    end
    
end

%% Export data

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