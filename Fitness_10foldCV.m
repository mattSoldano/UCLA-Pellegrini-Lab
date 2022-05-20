pause
clear; clc; close;
%%%%% Load ProteomicsPSLR workspace if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Name

savename1 = 'Proteome_PSLR_predict_vs_actual_traits_append_10foldCV.png';
savename2 = 'Proteome_PSLR_covariance_append_10foldCV.png';
exportName = 'Proteome_PSLR_predict_vs_actual_stats_append_10foldCV.xls';
savename3 = 'MainComps2traits_10foldCV.png';
savename4 = 'PSLR_Pvalues_append_Turk_10foldCV.png';
rescale = 0;
% _Rescaled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
methylT = readtable('MethMatrixCG-HumanBloodFitness96manual_20-100.csv', 'VariableNamingRule', 'preserve');
methylData = table2array(methylT(:, 2:end));
methylSites = methylT.Properties.VariableNames(2:end);
samplenames = methylT.sampleId;
%% Import traits and import the sample names we have methylation values for

T = readtable('Traits_ProsperProteome103_120821mjt.csv', 'VariableNamingRule', 'preserve');
Traits = T.Properties.VariableNames(2:end);
traitdata = table2array(T(:, 2:end));
traitnames = T.sampleID;

%% Rescale the data, so all traits range from 0 - 1 per sample

if rescale == 1
    for i = 1:width(traitdata)
        traitdata(:,i) = normalize(traitdata(:,i));
    end
end
%% Compare names from metadata and methylation value names to make sure they are the same

traitIdx = zeros(length(samplenames),1);

for kit=1:numel(samplenames)
    if length(find(traitnames == samplenames(kit))) == 1
        traitIdx(kit) = find(traitnames == samplenames(kit));
    end
end

%% Declare X and Y for PSLR
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

X = methylData;
Y = traitdata(traitIdx, :);
[mainA, mainB, mainU, mainV, mainBETA, mainPCTVAR, mainMSE, mainStats] = plsregress(X, Y, 20, 'CV', 10);
Ycat = horzcat(Y, mainV);

%% PSLR LOOCV
uniq = unique(samplenames);
numUniqueSamples = numel(uniq);
numsamples = height(Y);

for numpc = 1:20
    numpc
    pause(120)
    for k=1:numUniqueSamples
        idx = [1:numsamples]';
        id = table(idx, samplenames);
        clear removedID
        % removing one individual (both methylation and age values) logic for lasso function
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
        yTrainHorzcat = Ycat(id.idx, :);
        
        [a,b,u,v,beta,pctvar,mse,stats] = plsregress(XTrain,yTrainHorzcat, numpc, 'CV', 10);
        for dups = 1:numel(testSample)
            yfit(testSample(dups),:) = [ones(1,1) X(testSample(dups),:)] * beta;
        end
    end
    [Temp_modelcorr, Temp_pVal] = corr(Ycat,yfit);
    modelcorr(:,numpc) = diag(Temp_modelcorr);
    pVal(:,numpc) = diag(Temp_pVal);
end

%% Plot the lasso predictions compared to the actual Comp for a sanity check

MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
    ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];

%% predict traits

for t = 1:width(mainU)
    cornamesU(t,1) = sprintf("U%d", t);
    cornamesV(t,1) = sprintf("V%d", t);
    cornames(t,1) = sprintf("Comp. %d", t);
end
vbls_comps = [Traits,cornames'];
preditTraits1 = figure(1);
predictVSactual = corr(yfit, Ycat);
heatmap(vbls_comps, vbls_comps, predictVSactual, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Predicted traits 10foldCV');
ylabel('Traits');

for errorall = 1:width(Ycat)
    maeall(:,errorall) = mean(abs(Ycat(:,errorall) - yfit(:,errorall)));
end

%% show components vs traits

preditTraits2 = figure(2);
heatmap(cornames, vbls_comps, modelcorr, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Components');
ylabel('Traits and Appended Components 10foldCV');

%% Plot of correlations of components to traits

comps2traits = figure(3);
mainComp2traits = corr(mainV,Ycat);
heatmap(cornamesV, vbls_comps, mainComp2traits', 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Main Components (V)');
ylabel('Traits and Appended Main Components 10foldCV');

%% Save corr values and p values

pValHeatMap = figure(4);
heatmap(cornames, vbls_comps, log10(pVal), 'ColorLimits', [-5, 0], 'Colormap', MYmap);
xlabel('Predicted Components (log10(P Values))');
ylabel('Traits and Appended Main Components 10foldCV)');

%% Save corr values and p values

writecell({'Correlation Values'}, exportName, 'Range', 'A1');
writecell({cornames{:}}, exportName, 'Range', 'B1')
writecell({vbls_comps{:}}', exportName, 'Range', 'A2')
writematrix(modelcorr, exportName, 'Range', 'B2')
writecell({'P-Values'}, exportName, 'Sheet', 2, 'Range', 'A1')
writecell({cornames{:}}, exportName, 'Sheet', 2, 'Range', 'B1')
writecell({vbls_comps{:}}', exportName, 'Sheet', 2, 'Range', 'A2')
writematrix(pVal, exportName, 'Sheet', 2, 'Range', 'B2')
writecell({'MAE:'}, exportName, 'Sheet', 3, 'Range', 'B1')
writecell({vbls_comps{:}}, exportName, 'Sheet', 3, 'Range', 'A2')
writematrix(maeall, exportName, 'Sheet', 3, 'Range', 'B2')
%% Save png

saveas(preditTraits1, savename1);
saveas(preditTraits2, savename2);
saveas(comps2traits, savename3);
saveas(pValHeatMap, savename4);
save('ProteomicsPSLR_10foldCV');
