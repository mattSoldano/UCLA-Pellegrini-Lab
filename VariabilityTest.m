pause
clear; clc; close;
%% Export names
exportName = 'FitnessTomMarco_PCVariability.xlsx';

%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
methylT = readtable('MethMatrixCG-HumanBuccalFitness192_10-100.csv', 'VariableNamingRule', 'preserve');
methylData = table2array(methylT(:, 2:end));
methylSites = methylT.Properties.VariableNames(2:end);
samplenames = methylT.sampleID;
%% Import traits and import the sample names we have methylation values for

T = readtable('Traits_HumanBuccalFitness104x173_040422mjt.csv', 'VariableNamingRule', 'preserve');
TM = readtable('TraitsAndDates_MarcoTomBuccal88_042722mjt', 'VariableNamingRule', 'preserve');
TMage = TM.age;
marco = find(TMage == 33); tom = find(TMage ~= 33);
MarcoID = TM.sampleID(marco);
TomID = TM.sampleID(tom);
alltraitnames = [T.subjectID; MarcoID; TomID];

%% Compare names from metadata and methylation value names to make sure they are the same
traitIdx = zeros(length(samplenames),1);

for kit=1:numel(samplenames)
    if length(find(strcmp(alltraitnames,samplenames{kit}))) == 1
        allIdx(kit) = find(strcmp(alltraitnames,samplenames{kit}));
        if length(find(strcmp(alltraitnames,samplenames{kit}))) == 1 ...
                && isempty(find(strcmp(MarcoID,samplenames{kit}))) ...
                && isempty(find(strcmp(TomID,samplenames{kit})))
            traitIdx(kit) = find(strcmp(alltraitnames,samplenames{kit}));
        elseif length(find(strcmp(MarcoID,samplenames{kit}))) == 1
            marcoSid(kit) = find(strcmp(MarcoID,samplenames{kit}));
        elseif length(find(strcmp(TomID,samplenames{kit}))) == 1
            tomSid(kit) = find(strcmp(TomID,samplenames{kit}));
        end
    end
end
%% PCA

allId = find(allIdx ~= 0);
fitId = find(traitIdx ~= 0);
marcoSid = find(marcoSid ~= 0);
tomSid = find(tomSid ~= 0);
[coefficients, allScores] = pca(methylData(allId,:));

marcoScores = methylData(marcoSid,:) * coefficients;
tomScores = methylData(tomSid,:) * coefficients;
allScores = methylData(fitId, :) * coefficients;
%% plots

cornames(:,1) = 1:width(allScores);

load('SigPCs.mat')

for p = 1:numel(sigPCs)
    fig = figure(1);
    x = [allScores(:,sigPCs(p)); tomScores(:,sigPCs(p)); marcoScores(:,sigPCs(p))];
    l1 = repmat({'Fitness Study'}, height(allScores), 1);
    l2 = repmat({'Tom'}, height(tomScores), 1);
    l3 = repmat({'Marco'}, height(marcoScores), 1);
    l = [l1; l2; l3];
    violinplot(x, l);
    
    xl = sprintf("PC %.0f", sigPCs(p));
    xlabel(xl);
    ylabel("Component Value");
    savename = sprintf('boxplots_PC%d.png', sigPCs(p));
    
    saveas(fig, savename);
    close all;
    clear fig;
end