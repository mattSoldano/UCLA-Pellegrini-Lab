clear; clc; close;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Name
saveName = 'RawTraitData';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in the methylation matrix
methylT = readtable('MethMatrixCG-HumanBuccalTurk95_10-100.csv', 'VariableNamingRule', 'preserve');
methylCoords = methylT.Properties.VariableNames(2:end);
samplenames = methylT.sampleID;
methyldata = table2array(methylT(:,2:end));

%% Import traits and import the sample names we have methylation values for

T = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');
Vblnames = T.Properties.VariableNames(3:end);
traitnames = T.No;
% Same IDs as the trait metadata excel file
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
kitid = kitid';

%% Filter out names that we don't have methylation or metadata on
traitIdx = kitid(1:95);
% all valid indexes for finding the right PC1s
%% Declare X and Y for canoncorr
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

X = methyldata;

Y = data(traitIdx, :);

%% canoncorr

[A, B, r, U, V, stats] = canoncorr(X,Y);

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
