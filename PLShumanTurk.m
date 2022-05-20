
pause
clear; clc; close;
%%%%% Load PLShumanTurk workspace if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Name

savename1 = 'PSLR_methylComps_traitComps_Turk.png';
savename2 = 'PSLR_methylComps_traits_Turk.png';
% exportName = 'PSLR_predict_vs_actual_stats_Turk.xls';
savename3 = 'PSLR_traits_predictedTraits_Turk.png';
rescale = 0;
% _Rescaled
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

%% Declare X and Y for PSLR
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

X = methyldata;
Y = traitdata(traitIdx, :);
[mainA, mainB, mainU, mainV, mainBETA, mainPCTVAR, mainMSE, mainStats] = plsregress(X, Y, 20, 'CV', 10);
% Ycat = horzcat(Y, mainV);

%% PSLR LOOCV
uniq = unique(samplenames);
numUniqueSamples = numel(uniq);
numsamples = height(Y);

for numpc = 1:20
    numpc
    pause(15);
    for k=1:numUniqueSamples
        idx = [1:numsamples]';
        id = table(idx, samplenames);
        clear testSample
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
        YTrain = Y(id.idx, :);
%         yTrainHorzcat = Ycat(id.idx, :);
        
        [a,b,u,v,beta,pctvar,mse,stats] = plsregress(XTrain,YTrain, numpc);
        for dups = 1:numel(testSample)
            yfit(testSample(dups),:) = [ones(1,1) X(testSample(dups),:)] * beta;
        end
    end
    [Temp_modelcorr, Temp_pVal] = corr(Y,yfit);
    modelcorr(:,numpc) = diag(Temp_modelcorr);
    pVal(:,numpc) = diag(Temp_pVal);
end

%% Plot the lasso predictions compared to the actual Comp for a sanity check

MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
    ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];

%% predicted traits vs traits

for t = 1:width(mainU)
    cornamesU(t,1) = sprintf("U%d", t);
    cornamesV(t,1) = sprintf("V%d", t);
    cornames(t,1) = sprintf("Comp. %d", t);
end

mComp_tComp = figure(1);
mCtC = corr(mainV,mainU);
heatmap(cornamesU, cornamesV, mCtC, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Methyl Components (U)');
ylabel('Trait Components (V)');

for errorall = 1:width(Y)
    maeall(:,errorall) = mean(abs(Y(:,errorall) - yfit(:,errorall)));
end

%% show predicted components vs traits

mComp_T = figure(2);
mCt = corr(Y, mainU);
heatmap(cornamesU, Vblnames, mCt, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Methyl Components (U)');
ylabel('Traits');

%% Plot of correlations of components to traits

T_pT = figure(3);
predictVSactual = corr(Y, yfit);
heatmap(Vblnames, Vblnames, predictVSactual, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Predicted traits');
ylabel('Traits');

%% Save corr values and p values

% pValHeatMap = figure(4);
% heatmap(cornames, Vblnames, log10(pVal), 'ColorLimits', [-5, 0], 'Colormap', MYmap);
% xlabel('Predicted Components (log10(P Values))');
% ylabel('Traits and Appended Main Components (log10(P Values))');

%% P-Value, Corr and MAE table
% writecell({'Correlation Values'}, exportName, 'Range', 'A1');
% writecell({cornames{:}}, exportName, 'Range', 'B1')
% writecell({Vblnames{:}}', exportName, 'Range', 'A2')
% writematrix(modelcorr, exportName, 'Range', 'B2')
% writecell({'P-Values'}, exportName, 'Sheet', 2, 'Range', 'A1')
% writecell({cornames{:}}, exportName, 'Sheet', 2, 'Range', 'B1')
% writecell({Vblnames{:}}', exportName, 'Sheet', 2, 'Range', 'A2')
% writematrix(pVal, exportName, 'Sheet', 2, 'Range', 'B2')
% writecell({'MAE:'}, exportName, 'Sheet', 3, 'Range', 'A2')
% writecell({Vblnames{:}}, exportName, 'Sheet', 3, 'Range', 'B1')
% writematrix(maeall, exportName, 'Sheet', 3, 'Range', 'B2')
%% Save png

saveas(mComp_tComp, savename1);
saveas(mComp_T, savename2);
saveas(T_pT, savename3);
save('PLShumanTurk_noAppend');