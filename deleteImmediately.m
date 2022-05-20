
clear; clc; close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Name
exportName = 'CCmethylProbes_LOOCV.xls';
figNameV = 'CheckLOOCV_CCA_V.png';
figNameU = 'CheckLOOCV_CCA_U.png';
rescale = 0;
% _Rescaled
% _OnlyBestProbes
heatmapName = 'CC_Heatmap_LOOCV.png';
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
                if isequal(uniq{k},id.samplenames{sa})
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
        % I think we want to know how related Y is to U
        
        
        Vhat(testSample(:), :) = Y(testSample(:),:) * B;
        % I shouldn't be using their trait data to predict the components,
        % I need to predict the components without trait data
        
end

%% Plot the lasso predictions compared to the actual Comp for a sanity check

    modelcorrU = zeros(width(U),1);
    modelcorrV = zeros(width(V),1);
    fun = @(x) sprintf('%0.4f', x);
    % format the decimals of the correlation matrix

    for t = 1:width(traitdata)
        cornamesU(t,1) = sprintf("U%d", t);
        cornamesV(t,1) = sprintf("V%d", t);
        cornames(t,1) = sprintf("Comp. %d", t);
    end
    
    checkV = figure(1);
    
    for chV = 1:height(cornamesV)
        subplot(4,7,chV)
        scatter(mainV(:,chV),Vhat(:,chV));
        lablx = append('Actual ', cornamesV(chV));
        lably = append('Predicted ', cornamesV(chV));
        xlabel(lablx);
        ylabel(lably);
        modelcorrV(chV) = corr(mainV(:,chV), Vhat(:,chV));     
    end
    modelcorrV = cellfun(fun, num2cell(modelcorrV(:)), 'UniformOutput', 0);
    
    checkU = figure(2);
    
    for chU = 1:height(cornamesU)
        subplot(4,7,chU)
        scatter(mainU(:,chU),Uhat(:,chU));
        lablx = append('Actual ', cornamesU(chU));
        lably = append('Predicted ', cornamesU(chU));
        xlabel(lablx);
        ylabel(lably);
        modelcorrU(chU) = corr(mainU(:,chU), Uhat(:,chU));
    end
    modelcorrU = cellfun(fun, num2cell(modelcorrU), 'UniformOutput',0);

    %% Make the data ready for saving
    
    goodProbesIdxT = nan(numsamples, width(U));
    for PT = 1:width(U)
        samplename = uniq{PT};
        goodProbesIdxT(1:height(goodProbesIdx(:).(samplename)), PT) = goodProbesIdx(:).(samplename);
        
    end
    
    UniqueProbes = unique(goodProbesIdxT);
    UniqueProbes = rmmissing(UniqueProbes);
    
    probeCount = zeros(height(UniqueProbes), 1);
    for f = 1:height(UniqueProbes)
        probeCount(f, 1) = sum(goodProbesIdxT == UniqueProbes(f), 'all');
    end
    %% Save everything

    saveas(checkV, figNameV);
    saveas(checkU, figNameU);
    
    correlationsCell = [[{"Component"}; cellstr(cornamesU)], [{"R Value"}; modelcorrU(:)]...
        [{"Component"}; cellstr(cornamesV)], [{"R Value"}; modelcorrV(:)]];
    
    writecell(correlationsCell, exportName, 'Sheet', 1);
    
    
    writetable(array2table(goodProbesIdxT, 'VariableNames', cornames), exportName, 'Sheet', 2);
    writetable(array2table(UniqueProbes), exportName, 'Sheet', 3, 'Range', 'A1');
    writecell({'Frequency of Probe'}, exportName, 'Sheet', 3, 'Range', 'B1')
    writematrix(probeCount, exportName, 'Sheet', 3, 'Range', 'B2')
    writecell({'# Unique Probes'}, exportName, 'Sheet', 3, 'Range', 'C1')
    writematrix(numel(UniqueProbes), exportName, 'Sheet', 3, 'Range', 'C2')
    
%% %%%%%%%%%% Heatmap Section %%%%%%%%%%%%%%%%

CC_heatmap = figure(3);

MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
    ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];

modelcorr = corr(Y, Vhat);
heatmap(cornamesV, Vblnames, modelcorr, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Canonical Components');
ylabel('Traits');

saveas(CC_heatmap, heatmapName)

%% Export the probes with coefficients and the canonical correlation values

% writecell(goodProbes, 'CCmethylProbes.csv');

% writetable(array2table(V), 'CcScoresForTraits.csv');