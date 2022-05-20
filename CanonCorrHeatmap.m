clear; clc; close;

%%%%%%% UNIX Commands for Files too large for readtable%%%%%%%%%%%%%%%%%%%
%{
tail -n +2 filename > newfile
cut -f 2- -d "," newfile > finalfile
%^ this is to remove the headers and the column with sample names

head -n 1 filename > newProbenameFile
cut -f 1 filename > newSamplenameFile
%^ These are to extract the headers into a new file

then in Matlab use

dlmread('finalfile',',')
methylCoords = string(readcell('MethMatrixProbeNames.csv'));

Unix example for this script:

tail -n +2 MethMatrixCG-HumanBuccalTurk95_10-100_Full.csv > MethMatrixTurkNoHeaders.csv
cut -f 2- -d "," MethMatrixTurkNoHeaders.csv > MethMatrixTurk.csv

head -n 1 MethMatrixCG-HumanBuccalTurk95_10-100_Full.csv > MethMatrixProbeNames.csv
cut -f 1 -d "," MethMatrixCG-HumanBuccalTurk95_10-100_Full.csv > MethMatrixSampleNames.csv
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Name
saveName = 'RawTraitData';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
noSortData = dlmread('MethMatrixTurk.csv', ',');
methylCoords = string(readcell('MethMatrixProbeNames.csv'));
methylCoords = methylCoords(2:end);
originalnames = readtable('MethMatrixSampleNames.csv');

%% Rename samplenames
turk = 'Turk-';
r = 'R';
fixnames = erase(originalnames.sampleID, turk);
fixnames = erase(fixnames, r);
for s = 1:height(fixnames)

    if strlength(fixnames(s)) == 1
        fixnames(s) = append('00', fixnames(s));
    elseif strlength(fixnames(s)) == 2
        fixnames(s) = append('0', fixnames(s));
    end
    
    if contains(fixnames(s), 'C') == 0
       fixnames(s) = append('E', fixnames(s));
    end

end

samplenames = sort(fixnames);

%% Unique for Turksample sort methylation values by samplenames

id = 1;
samplenameIDX = zeros(length(samplenames),1);
while id <= numel(samplenames)
    if length(find(strcmp(fixnames,samplenames{id}))) == 1
        samplenameIDX(id) = find(strcmp(fixnames,samplenames{id}));
        id = id + 1;
    elseif length(find(strcmp(fixnames,samplenames{id}))) == 2
        temp = find(strcmp(fixnames,samplenames{id}));
        samplenameIDX(id) = temp(1);
        id = id + 1;
        samplenameIDX(id) = temp(2);
        id = id + 1;
    end
end

% samplenameIDX = zeros(length(samplenames),1);
% for id = 1:numel(samplenames)
%         samplenameIDX(id) = find(strcmp(fixnames,samplenames{id}));
% end

methyldata = noSortData(samplenameIDX, :);

%% Import traits and import the sample names we have methylation values for

T = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');
Vblnames = T.Properties.VariableNames(3:end);
traitnames = T.No;
data = knnimpute(table2array(T(:, 3:end)));
desiredTraits = [1 2 7:28 30 31];
data = data(:, desiredTraits);
Vblnames = Vblnames(desiredTraits);

% %% Rescale the data, so all traits range from 0 - 1 per sample
% 
% rescaleData = zeros(height(data), width(data));
% 
% for i = 1:width(data)
% data(:,i) = rescale(data(:,i));
% end

%% Compare names from metadata and methylation value names to make sure they are the same

kitid = zeros(length(traitnames),1);

for kit=1:numel(samplenames)
    if length(find(strcmp(traitnames,samplenames{kit}))) == 1
        kitid(kit) = find(strcmp(traitnames,samplenames{kit}));
    end
end

%% Filter out names that we don't have methylation or metadata on
traitIdx = kitid(1:numel(samplenames));
% all valid indexes for finding the right PC1s
%% Declare X and Y for canoncorr
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

X1 = methyldata;

Y1 = data(traitIdx, :);

%% canoncorr

[A, B, r, U, V, stats] = canoncorr(X1,Y1);
% V is what we are doing the analysis on (the canonical components)

%% Plot CannonCorr

correlation = figure(1);

for k = 1:width(data)
cornamesX(k,1) = sprintf("U%d", k);
cornamesY(k,1) = sprintf("V%d", k);
subplot(8, 4, k);
plot(U(:,k), V(:,k), '.');
xlabel(cornamesX(k));
ylabel(cornamesY(k));
end

pngSave = append(saveName, '_CorrPlot_Raw_LimitedTraits_noSort.png');
saveas(correlation, pngSave);

% Export the components for each sample using the V matrix

tableV = array2table(V, 'VariableNames', cornamesY);
componentSave = append(saveName, '_ccorrComponents_Raw_LimitedTraits_noSort.csv');
writetable(tableV, componentSave);

%%%%%%%%%%%% Heatmap Section %%%%%%%%%%%%%%%%

%% Determining Correlation
% Make sure X has the correct amount of samples and that the samples are
% properly aligned
  
        Y = data(traitIdx,:);
        X = V;
        modelcorr = corr(X,Y);

CC_heatmap = figure(2);

heatmap(cornamesY, Vblnames, modelcorr, 'ColorLimits', [-1, 1], 'Colormap', turbo);
xlabel('Canonical Components');
ylabel('Traits');

saveas(CC_heatmap, 'CC_Heatmap_Raw_LimitedTraits_noSort.png')