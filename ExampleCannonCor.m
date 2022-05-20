clear; clc; close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Name
saveName = 'Example';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in the data (carbig file) in replacement of the methylation and trait matrices
load carbig
data = [Displacement Horsepower Weight Acceleration MPG];
nans = sum(isnan(data),2) > 0;

%% Compare names from metadata and methylation value names to make sure they are the same
% Not needed with this dataset 

% kitid = zeros(length(traitnames),1);
% 
% for kit=1:numel(samplenames)
%     if length(find(strcmp(traitnames,samplenames{kit}))) == 1
%         kitid(kit) = find(strcmp(traitnames,samplenames{kit}));
%     end
% end
% kitid = kitid';

%% Filter out names that we don't have methylation or metadata on

% Not needed

% traitIdx = kitid(1:95);
% all valid indexes for finding the right PC1s
%% Declare X and Y for canoncorr
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

X = methyldata;

Y = data(traitIdx, :);

%% canoncorr

[A, B, r, U, V, stats] = canoncorr(X,Y);
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
componentSave = append(saveName, '_ccorrComponents.csv');
writetable(tableV, componentSave);

%%%%%%%%%%%% Heatmap Section %%%%%%%%%%%%%%%%

%% Determining Correlation
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

% modelcorr = zeros(width(Vblnames), width(cornamesY));
% 
% for cc = 1:width(V)
%     for t = 1:width(data)
        
        Y = data(traitIdx,:);
        X = V;
        
        

        modelcorr = corr(X,Y);

%     end
% end

CC_heatmap = figure(2);

heatmap(cornamesY, Vblnames, modelcorr, 'ColorLimits', [-1, 1], 'Colormap', turbo);
xlabel('Canonical Components');
ylabel('Traits');

saveas(CC_heatmap, 'CC_Heatmap_Example.png')