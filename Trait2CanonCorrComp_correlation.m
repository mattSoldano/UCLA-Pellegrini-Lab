clear; clc; close;

%% Import Canonical Components
CCtable = readtable('RawTraitData_ccorrComponents_Raw_LimitedTraits_noSort.csv', 'VariableNamingRule', 'preserve');
CCnames = CCtable.Properties.VariableNames;
CCs = table2array(CCtable);

%% Import identifiers for canonical correlation/methylation indices
methylT = readtable('MethMatrixCG-HumanBuccalTurk95_10-100.csv', 'VariableNamingRule', 'preserve');
samplenames = methylT.sampleID;

%% Import the traits and the identifiers for the trait indices
T = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');
Vblnames = T.Properties.VariableNames(3:end);
traitnames = T.No;
% Same IDs as the trait metadata excel file
traitdata = knnimpute(table2array(T(:, 3:end)));
desiredTraits = [1 2 7:28 30 31];
traitdata = traitdata(:, desiredTraits);
Tvbls = Vblnames(desiredTraits);

%% Compare names from metadata and methylation value names to make sure they are the same

kitid = zeros(length(traitnames),1);

for kit=1:numel(samplenames)
    if length(find(strcmp(traitnames,samplenames{kit}))) == 1
        kitid(kit) = find(strcmp(traitnames,samplenames{kit}));
    end
end
kitid = kitid';

%% Filter out names that we don't have data on
traitIdx = kitid(1:numel(samplenames));

%% Declare X and Y for Correlation
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

modelcorr = zeros(width(Tvbls), width(CCnames));

for cc = 1:width(CCs)
    for t = 1:width(traitdata)
        
        Y = traitdata(traitIdx,t);
        X = CCs(:, cc);
        
        %% Determining Correlation

        modelcorr(t, cc) = corr(X,Y);

    end
end

CC_heatmap = figure(1);

heatmap(CCnames, Tvbls, modelcorr, 'ColorLimits', [-1, 1], 'Colormap', turbo);
xlabel('Canonical Components');
ylabel('Traits');

saveas(CC_heatmap, 'CC_Heatmap_Raw_LimitedTraits_noSort.png')