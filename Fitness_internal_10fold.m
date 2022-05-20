pause
clear; clc; close;
%%%%% Load ProteomicsPSLR workspace if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Name

savename1 = 'PSLR_methylComps_traitComps_Turk_internal_10foldCV_Neiman.png';
savename2 = 'PSLR_methylComps_traits_Turk_internal_10foldCV_Neiman.png';
savename3 = 'PSLR_traits_predictedTraits_Turk_internal_10foldCV_Neiman.png';
% savename1 = 'Proteome_PSLR_predict_vs_actual_traits_internal_10foldCV_Neiman.png';
% savename2 = 'Proteome_PSLR_covariance_internal_10foldCV_Neiman.png';
% exportName = 'Proteome_PSLR_predict_vs_actual_stats_internal_10foldCV_Neiman.xls';
% savename3 = 'MainComps2traits_internal_10foldCV_Neiman.png';
% savename4 = 'PSLR_Pvalues_Turk_internal_10foldCV_Neiman.png';
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
[mainA, mainB, mainU, mainV, mainBETA, mainPCTVAR, mainMSE, mainStats] = plsregress(X, Y, 20);
% Ycat = horzcat(Y, mainV);

%% PSLR LOOCV
numsamples = height(Y);
for numpc = 1:20
    numpc
    [a,b,u,v,beta,pctvar,mse,stats] = plsregress(X,Y, numpc, 'CV', numsamples - 1);
    yfit = [ones(numsamples,1) X] * beta;
    
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
xlabel('Methyl Components (U) internal 10foldCV');
ylabel('Trait Components (V) internal 10foldCV');

for errorall = 1:width(Y)
    maeall(:,errorall) = mean(abs(Y(:,errorall) - yfit(:,errorall)));
end

%% show predicted components vs traits

mComp_T = figure(2);
mCt = corr(Y, mainU);
heatmap(cornamesU, Traits, mCt, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Methyl Components (U) internal 10foldCV');
ylabel('Traits');

%% Plot of correlations of components to traits

T_pT = figure(3);
predictVSactual = corr(Y, yfit);
heatmap(Traits, Traits, predictVSactual, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Predicted traits internal 10foldCV');
ylabel('Traits');

%% Save png

saveas(mComp_tComp, savename1);
saveas(mComp_T, savename2);
saveas(T_pT, savename3);

% %% Plot the lasso predictions compared to the actual Comp for a sanity check
% 
% MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
%     ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];
% 
% %% predict traits
% 
% for t = 1:width(mainU)
%     cornamesU(t,1) = sprintf("U%d", t);
%     cornamesV(t,1) = sprintf("V%d", t);
%     cornames(t,1) = sprintf("Comp. %d", t);
% end
% 
% % vbls_comps = [Traits,cornames'];
% preditTraits1 = figure(1);
% predictVSactual = corr(yfit, Y);
% heatmap(Traits, Traits, predictVSactual, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
% xlabel('Predicted traits internal 10foldCV');
% ylabel('Traits');
% 
% for errorall = 1:width(Y)
%     maeall(:,errorall) = mean(abs(Y(:,errorall) - yfit(:,errorall)));
% end
% 
% %% show components vs traits
% 
% preditTraits2 = figure(2);
% heatmap(cornames, Traits, modelcorr, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
% xlabel('Components');
% ylabel('Traits internal 10foldCV');
% 
% %% Plot of correlations of components to traits
% 
% comps2traits = figure(3);
% mainComp2traits = corr(mainV,Y);
% heatmap(cornamesV, Traits, mainComp2traits', 'ColorLimits', [-1, 1], 'Colormap', MYmap);
% xlabel('Main Components (V)');
% ylabel('Traits internal 10foldCV');
% 
% %% Save corr values and p values
% 
% pValHeatMap = figure(4);
% heatmap(cornames, Traits, log10(pVal), 'ColorLimits', [-5, 0], 'Colormap', MYmap);
% xlabel('Predicted Components (log10(P Values))');
% ylabel('Traits internal 10foldCV)');
% 
% %% Save corr values and p values
% 
% writecell({'Correlation Values'}, exportName, 'Range', 'A1');
% writecell({cornames{:}}, exportName, 'Range', 'B1')
% writecell({Traits{:}}', exportName, 'Range', 'A2')
% writematrix(modelcorr, exportName, 'Range', 'B2')
% writecell({'P-Values'}, exportName, 'Sheet', 2, 'Range', 'A1')
% writecell({cornames{:}}, exportName, 'Sheet', 2, 'Range', 'B1')
% writecell({Traits{:}}', exportName, 'Sheet', 2, 'Range', 'A2')
% writematrix(pVal, exportName, 'Sheet', 2, 'Range', 'B2')
% writecell({'MAE:'}, exportName, 'Sheet', 3, 'Range', 'B1')
% writecell({Traits{:}}, exportName, 'Sheet', 3, 'Range', 'A2')
% writematrix(maeall, exportName, 'Sheet', 3, 'Range', 'B2')
% %% Save png
% 
% saveas(preditTraits1, savename1);
% saveas(preditTraits2, savename2);
% saveas(comps2traits, savename3);
% saveas(pValHeatMap, savename4);
% save('ProteomicsPSLR_internal_10foldCV');
