
clear; clc; close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Name
exportName = 'CCmethylProbes_methyl_LOOCV.xls';
figNameALL = 'Methyl_CCA_Heatmaps.png';
rescale = 0;
% _Rescaled
% _OnlyBestProbes
% heatmapNameLOO = 'Check_methyl_LOOCV_CCA_UvsV.png';
% heatmapNameCCA = 'Check_methyl_CCA_UvsV.png';
% heatmapNameLOO_Traits = 'CC_traits_methyl_Heatmap_LOOCV.png';
% heatmapNameCCA_Traits = 'CC_traits_methyl_Heatmap_CCA.png';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
methyldata = dlmread('MethMatrixTurk.csv', ',');
methylCoords = readcell('MethMatrixProbeNames.csv');
methylCoords = methylCoords(2:end);
samplenames = readcell('MethMatrixSampleNames.csv');
samplenames = (erase(samplenames(2:end), 'R'));
samplenames = regexprep(samplenames, '-', '_');
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
% change the gender value --> Female = 1, male = 0
T.Gender = cell2mat(T.Gender);
desiredTraits = [3 4 8:29 31 32];
traitnames = readcell('MethMatrixTraitNames.csv');
traitnames = traitnames(2:end);
traitdata = table2array(T(:, desiredTraits));
Vblnames = Vblnames(desiredTraits);

%% Remove sample C40 because of NaN values
C40t = find(strcmp(traitnames, "C40"));
traitnames(C40t) = [];
traitdata(C40t,:) = [];

C40m = find(strcmp(samplenames, "Turk_C40"));
methyldata(C40m, :) = [];
samplenames(C40m, :) = [];

%% Rename traitnames
turk = 'Turk_';
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

%% Main CCA

[mainA, mainB, mainR, mainU, mainV] = canoncorr(X, Y);
% find best probes
mainProbes = find(mainA(:,1) ~= 0);

%% canoncorr LOOCV
uniq = unique(samplenames);
numUniqueSamples = numel(uniq);
numsamples = height(Y);
Uhat = zeros(height(mainU),width(mainU));
Vhat = zeros(height(mainV),width(mainV));

for k=1:numUniqueSamples
    idx = [1:numsamples]';
    id = table(idx, samplenames);
    clear removedID
    % removing one individual (both methylation and age values) logic for lasso function
    for j = 1:numsamples
        if isequal(uniq{k},id.samplenames{j}) == true
            testSample(1) = j;
            id(j,:) = [];
            for sa = 1:height(id.samplenames)
                if isequal(uniq{k},id.samplenames{sa}) ... 
%                         && j ~= sa
%                     testSample(2) = sa;
                    testSample(2) = sa + 1;
                    id(sa,:) = [];
                    break
                end
            end
            break
        else
        end
    end
    
    XTrain = X(id.idx, :);
    yTrain = Y(id.idx, :);
    
    [A, B, r, U, V, stats] = canoncorr(XTrain,yTrain);
    
    %% Use the cooefficients from the training data to predict the sample that was left out
    samplename = uniq{k};
    goodProbesIdx(:).(samplename) = find(A(:,1) ~= 0);
    
    Uhat(testSample(:), :) = X(testSample(:), :) * A;
    Vhat(testSample(:), :) = Y(testSample(:),:) * B;
    
end

%% Plot the lasso predictions compared to the actual Comp for a sanity check

fun = @(x) sprintf('%0.4f', x);
% format the decimals of the correlation matrix

for t = 1:width(traitdata)
    cornamesU(t,1) = sprintf("U%d", t);
    cornamesV(t,1) = sprintf("V%d", t);
    cornames(t,1) = sprintf("Comp. %d", t);
end
% make the arrays with the component names 

% %% Scatterplot
% heatmap colors

% scat = figure(1);
% for sca = 1:width(Uhat)
%     subplot(6,5,sca)
%     scatter(Vhat(:,sca),Uhat(:,sca));
%     Vlabel = cornamesV(sca);
%     Ulabel = cornamesU(sca);
%     xlabel(Vlabel);
%     ylabel(Ulabel);
% end

MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
    ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];

check = figure(1);
% checkV = figure(1);
subplot(2,2,1)
modelcorrLOO = corr(Uhat,Vhat);
heatmap(cornamesV, cornamesU, modelcorrLOO, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Predicted V (Vhat)');
ylabel('Predicted U (Uhat)');

% checkCCA = figure(2);
subplot(2,2,2)
modelcorrCCA = corr(mainU, mainV);
heatmap(cornamesV, cornamesU, modelcorrCCA, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('mainV');
ylabel('mainU');

%% Make the data ready for saving

modelcorrCCA_E = zeros(height(modelcorrCCA),1);
modelcorrE = zeros(height(modelcorrLOO),1);

for h = 1:height(modelcorrLOO)
    modelcorrE(h,1) = modelcorrLOO(h,h);
end
for h = 1:height(modelcorrCCA)
    modelcorrCCA_E(h,1) = modelcorrCCA(h,h);
end

modelcorrExportCCA = cellfun(fun, num2cell(modelcorrCCA_E), 'UniformOutput',0);
modelcorrExportLOO = cellfun(fun, num2cell(modelcorrE), 'UniformOutput',0);

goodProbesIdxT = nan(numsamples, numUniqueSamples);
for PT = 1:numUniqueSamples
    samplename = uniq{PT};
    goodProbesIdxT(1:height(goodProbesIdx(:).(samplename)), PT) = goodProbesIdx(:).(samplename);
end

UniqueProbes = unique(goodProbesIdxT);
UniqueProbes = rmmissing(UniqueProbes);

probeCount = zeros(height(UniqueProbes), 1);
mup = 1;

for f = 1:height(UniqueProbes)
    probeCount(f, 1) = sum(goodProbesIdxT == UniqueProbes(f), 'all');
    if probeCount(f,1) >= 70
        MostUsedProbesIDX(mup,1) = f;
        mup = mup + 1;
    end
end

correlationsCell = [[{"Predicted"}; cellstr(cornames)], [{"R Value"}; modelcorrExportLOO(:)]...
    [{"Pure CCA"}; cellstr(cornames)], [{"R Value"}; modelcorrExportCCA(:)]];
%% %%%%%%%%%% Trait Heatmap Section %%%%%%%%%%%%%%%%

% CC_heatmapLOO = figure(3);
subplot(2,2,3)
modelcorrLOO = corr(Y, Vhat);
heatmap(cornamesV, Vblnames, modelcorrLOO, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Predicted Canonical Components');
ylabel('Traits');

% CC_heatmapCCA = figure(4);
subplot(2,2,4)
modelcorrCCA = corr(Y, mainV);
heatmap(cornamesV, Vblnames, modelcorrCCA, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Main Canonical Components');
ylabel('Traits');

%% Save everything

% PNGs
saveas(check, figNameALL);
% saveas(CC_heatmapLOO, heatmapNameLOO_Traits)
% saveas(CC_heatmapCCA, heatmapNameCCA_Traits)
% saveas(checkV, heatmapNameLOO);
% saveas(checkCCA, heatmapNameCCA);

% Spreadsheet
writecell(correlationsCell, exportName, 'Sheet', 1);
writetable(array2table(goodProbesIdxT, 'VariableNames', uniq), exportName, 'Sheet', 2);
writetable(array2table(UniqueProbes), exportName, 'Sheet', 3, 'Range', 'A1');
writecell({'Frequency of Probe'}, exportName, 'Sheet', 3, 'Range', 'B1')
writematrix(probeCount, exportName, 'Sheet', 3, 'Range', 'B2')
writecell({'# Unique Probes'}, exportName, 'Sheet', 3, 'Range', 'C1')
writematrix(numel(UniqueProbes), exportName, 'Sheet', 3, 'Range', 'C2')
writecell({'Most Used Probes'}, exportName, 'Sheet', 3, 'Range', 'D1')
writematrix(MostUsedProbesIDX, exportName, 'Sheet', 3, 'Range', 'D2')
