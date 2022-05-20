
clear; clc; close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Name
saveName = 'RawTraitData';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
methyldata = dlmread('MethMatrixTurk.csv', ',');
methylCoords = readcell('MethMatrixProbeNames.csv');
methylCoords = methylCoords(2:end);
samplenames = readcell('MethMatrixSampleNames.csv');
samplenames =(erase(samplenames(2:end), 'R'));
rescale = 0;
% _Rescaled
%% Import traits and import the sample names we have methylation values for

T = readtable('Turkish_Buccal_Traits2.csv', 'VariableNamingRule', 'preserve');
Vblnames = T.Properties.VariableNames;

for g = 1:height(T.Gender)
    if isequal(T.Gender(g), "M") == true
        T.Gender{g} = 0;
    else
        T.Gender{g} = 1;
    end
end
T.Gender = cell2mat(T.Gender);
desiredTraits = [3 4 8:29 31 32];
traitnames = readcell('MethMatrixTraitNames.csv');
traitnames = traitnames(2:end);
% traitdata = knnimpute(table2array(T(:, desiredTraits)));
% traitdata = table2array(T(:, desiredTraits));
Vblnames = Vblnames(desiredTraits);

%% Remove sample C40 because of NaN values
C40t = find(strcmp(traitnames, "C40"));
traitnames(C40t) = [];
traitdata = table2array(T(:, desiredTraits));
traitdata(C40t,:) = [];

C40m = find(strcmp(samplenames, "Turk-C40"));
methyldata(C40m, :) = [];
samplenames(C40m, :) = [];

%% Rename traitnames
turk = 'Turk-';
for s = 1:height(traitnames)
    traitnames{s} = append(turk, num2str(traitnames{s}));
end

%% Rescale the data, so all traits range from 0 - 1 per sample

if rescale == 1
    for i = 1:width(traitdata)
        traitdata(:,i) = normalize(traitdata(:,i));
    end
end

%% Compare names from metadata and methylation value names to make sure they are the same

traitIdx = zeros(length(samplenames),1);

for kit=1:numel(samplenames)
    if length(find(strcmp(traitnames,samplenames{kit}))) == 1
        traitIdx(kit) = find(strcmp(traitnames,samplenames{kit}));
    end
end

%% Declare X and Y for canoncorr
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

X = methyldata;

Y = traitdata(traitIdx, :);

%% canoncorr

[A, B, r, U, V, stats] = canoncorr(X,Y);
% V is what we are doing the analysis on (the canonical components)

%% Plot CannonCorr

correlation = figure(1);

for t = 1:width(traitdata)
    cornamesU(t,1) = sprintf("U%d", t);
    cornamesV(t,1) = sprintf("V%d", t);
    cornames(t,1) = sprintf("Comp. %d", t);
end
    
corrUV = corr(U,V);
MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
    ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];

heatmap(cornamesV, cornamesU, corrUV, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('V');
ylabel('U');

%%%%%%%%%%%% Heatmap Section %%%%%%%%%%%%%%%%

%% Determining Correlation
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

modelcorr = corr(Y,V);
% in the corr function, the first matrix corresponds to the rows (y axis) and
% the second matrix corresponds to the columns (x axis)

CC_heatmap = figure(5);

heatmap(cornamesV, Vblnames, modelcorr, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Canonical Components');
ylabel('Traits');

saveas(CC_heatmap, 'CC_Heatmap_Raw_LimitedTraits_noSort.png')

%% Export the probes with coefficients and the canonical correlation values

goodProbesIdx = find(A(:,1) ~= 0);
goodProbes = methylCoords(goodProbesIdx);
writecell(goodProbes, 'CCmethylProbes.csv');

writetable(array2table(V), 'CcScoresForTraits.csv');