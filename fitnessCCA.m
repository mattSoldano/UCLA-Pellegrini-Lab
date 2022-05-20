clear; clc; close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Name
exportName = 'CCproteomeProbes_LOOCV_Rescaled.xls';
% heatmapNameALL = 'Proteome_CCA_Heatmaps.png';
rescale = 1;
% _Rescaled
% _OnlyBestProbes
heatmapNameLOO = 'Checkproteome_heatmap_LOOCV_Rescaled.png';
heatmapNameCCA = 'Checkproteome_heatmap_CCA_Rescaled.png';
heatmapNameLOO_traits = 'CC_trait_proteome_Heatmap_LOOCV_Rescaled.png';
heatmapNameCCA_traits = 'CC_trait_proteome_Heatmap_CCA_Rescaled.png';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in the proteome matrices --> rename sampleIDs (manually in a txt file unsorted)
proteomeT = readtable('ProtMatrix103_120821mjt.csv', 'VariableNamingRule', 'preserve');
proteomeData = table2array(proteomeT(:, 2:end));
proteomeProteins = proteomeT.Properties.VariableNames(2:end);
samplenames = proteomeT.sampleID;

%% Import traits and import the sample names we have proteomeation values for

T = readtable('Traits_ProsperProteome103_120821mjt.csv', 'VariableNamingRule', 'preserve');
Traits = T.Properties.VariableNames(2:end);
traitdata = table2array(T(:, 2:end));
traitnames = T.sampleID;
% traitdata = knnimpute(table2array(T(:, desiredTraits)));

%% Rescale the data, so all traits range from 0 - 1 per sample

if rescale == 1
    if rescale == 1
        for i = 1:width(traitdata)
            traitdata(:,i) = normalize(traitdata(:,i));
        end
    end
end
%% Compare names from metadata and proteomeation value names to make sure they are the same

traitIdx = zeros(length(samplenames),1);

for kit=1:numel(samplenames)
    if length(find(traitnames == samplenames(kit))) == 1
        traitIdx(kit) = find(traitnames == samplenames(kit));
    end
end

%% Declare X and Y for canoncorr
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

X = proteomeData;
Y = traitdata(traitIdx, :);

%% Main CCA

[mainA, mainB, mainR, mainU, mainV] = canoncorr(X, Y);
% find best probes
mainProteins = find(mainA(:,1) ~= 0);

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
    % removing one individual (both proteomeation and age values) logic for lasso function
    for j = 1:numsamples
        if isequal(uniq(k),id.samplenames(j)) == true
            testSample(1) = j;
            id(j,:) = [];
            for sa = 1:height(id.samplenames)
                if isequal(uniq(k),id.samplenames(sa))
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
    samplename(k) = sprintf("Sample_%d", uniq(k));
    
    goodProteinsIdx(:).(samplename(k)) = find(A(:,1) ~= 0);
    
    Uhat(testSample(:), :) = X(testSample(:), :) * A;
    % I think we want to know how related Y is to U
    
    
    Vhat(testSample(:), :) = Y(testSample(:),:) * B;
    % I shouldn't be using their trait data to predict the components,
    % I need to predict the components without trait data
    
end

%% Plot the lasso predictions compared to the actual Comp for a sanity check

modelcorrExport = zeros(width(U),1);
modelcorrV = zeros(width(V),1);
fun = @(x) sprintf('%0.4f', x);
% format the decimals of the correlation matrix

for t = 1:width(U)
    cornamesU(t,1) = sprintf("U%d", t);
    cornamesV(t,1) = sprintf("V%d", t);
    cornames(t,1) = sprintf("Comp. %d", t);
end

% check = figure(1);
checkV = figure(1);
% subplot(2,2,1)
modelcorrLOO = corr(Uhat,Vhat);

MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
    ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];

heatmap(cornamesV, cornamesU, modelcorrLOO, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Predicted V (Vhat) Rescaled');
ylabel('Predicted U (Uhat) Rescaled');

checkCCA = figure(2);
% subplot(2,2,2)
modelcorrCCA = corr(mainU, mainV);
heatmap(cornamesV, cornamesU, modelcorrCCA, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('mainV Rescaled');
ylabel('mainU Rescaled');

%% Make the data ready for saving

modelcorrCCA_E = zeros(height(modelcorrCCA),1);
modelcorrLOO_E = zeros(height(modelcorrLOO),1);
for h = 1:height(modelcorrCCA)
    modelcorrCCA_E(h,1) = modelcorrCCA(h,h);
end
for h = 1:width(modelcorrLOO)
    modelcorrLOO_E(h,1) = max(modelcorrLOO(:,h));
end
modelcorrExportCCA = cellfun(fun, num2cell(modelcorrCCA_E), 'UniformOutput',0);
modelcorrExportLOO = cellfun(fun, num2cell(modelcorrLOO_E), 'UniformOutput',0);

goodProteinsIdxT = nan(numsamples, numUniqueSamples);
for PT = 1:numUniqueSamples
    goodProteinsIdxT(1:height(goodProteinsIdx(:).(samplename(PT))), PT) = goodProteinsIdx(:).(samplename(PT));
end

UniqueProteins = unique(goodProteinsIdxT);
UniqueProteins = rmmissing(UniqueProteins);

proteinCount = zeros(height(UniqueProteins), 1);
mup = 1;

for f = 1:height(UniqueProteins)
    proteinCount(f, 1) = sum(goodProteinsIdxT == UniqueProteins(f), 'all');
    if proteinCount(f,1) >= 70
        MostUsedProteinsIDX(mup,1) = f;
        mup = mup + 1;
    end
end

correlationsCell = [[{"Predicted Rescaled"}; cellstr(cornames)], [{"R Value"}; modelcorrExportLOO(:)]...
    [{"Pure CCA Rescaled"}; cellstr(cornames)], [{"R Value"}; modelcorrExportCCA(:)]];

%% %%%%%%%%%% Trait Heatmap Section %%%%%%%%%%%%%%%%

CC_heatmapLOO = figure(3);
% subplot(2,2,3)

modelcorrLOO = corr(Y, Vhat);
heatmap(cornamesV, Traits, modelcorrLOO, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Predicted Canonical Components Rescaled');
ylabel('Traits Rescaled');

CC_heatmapCCA = figure(4);
% subplot(2,2,4)

modelcorrCCA = corr(Y, mainV);
heatmap(cornamesV, Traits, modelcorrCCA, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Main Canonical Components Rescaled');
ylabel('Traits Rescaled');
%% Save everything

% saveas(check, heatmapNameALL);
saveas(checkV, heatmapNameLOO);
saveas(checkCCA, heatmapNameCCA);
saveas(CC_heatmapLOO, heatmapNameLOO_traits)
saveas(CC_heatmapCCA, heatmapNameCCA_traits)

writecell(correlationsCell, exportName, 'Sheet', 1);
writetable(array2table(goodProteinsIdxT, 'VariableNames', samplename), exportName, 'Sheet', 2);
writetable(array2table(UniqueProteins), exportName, 'Sheet', 3, 'Range', 'A1');
writecell({'Frequency of Protein Rescaled'}, exportName, 'Sheet', 3, 'Range', 'B1')
writematrix(proteinCount, exportName, 'Sheet', 3, 'Range', 'B2')
writecell({'# Unique Proteins Rescaled'}, exportName, 'Sheet', 3, 'Range', 'C1')
writematrix(numel(UniqueProteins), exportName, 'Sheet', 3, 'Range', 'C2')
writecell({'Most Used Proteins Rescaled'}, exportName, 'Sheet', 3, 'Range', 'D1')
writematrix(MostUsedProteinsIDX, exportName, 'Sheet', 3, 'Range', 'D2')