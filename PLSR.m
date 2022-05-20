
clear; clc; close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Name
exportName = 'CCmethylProbes_methyl_LOOCV_PSLR_PSLRappend.xls';
figNameALL = 'Methyl_PSLR_Figs_CompAppend.png';
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

%% Declare X and Y for PSLR
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

X = methyldata;
Y = traitdata(traitIdx, :);

%% Main PLSR --> XL = A, YL = B, XS = U, YS = V
[mainA, mainB, mainU, mainV, mainBETA, mainPCTVAR, mainMSE, mainStats] = plsregress(X, Y, 26, 'CV', 10);
% 
% %% Analyze PLSR
% W0 = mainStats.W ./ sqrt(sum(mainStats.W.^2,1));
% p = size(mainA,1);
% sumSq = sum(mainU.^2,1).*sum(mainB.^2,1);
% vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
% indVIPs = find(vipScore >= 1);
% find best probes

%% PSLR LOOCV
uniq = unique(samplenames);
numUniqueSamples = numel(uniq);
numsamples = height(Y);
Uhat = zeros(height(mainU),width(mainU*2));
Vhat = zeros(height(mainV),width(mainV*2));

% for t = 1:width(mainU)
%     cornamesU(t,1) = sprintf("U%d", t);
%     cornamesV(t,1) = sprintf("V%d", t);
%     cornames(t,1) = sprintf("Comp. %d", t);
% end

mainY = Y;
% Y = horzcat(Y, mainU);
yfit = zeros(height(mainY), width(mainY));
for numpc = 1:20
    numpc
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
        yTrain = mainY(id.idx, :);
        %     yTrainHorzcat = Y(id.idx, :);
        
        %     [A,B,U,V,BETA,PCTVAR,MSE,STATS] = plsregress(XTrain,yTrainHorzcat, 26, 'CV', 10);
        [a,b,u,v,beta,pctvar,mse,stats] = plsregress(XTrain,yTrain, numpc, 'CV', 10);
        for dups = 1:numel(testSample)
            yfit(testSample(dups),:) = [ones(1,1) X(testSample(dups),:)] * beta;
        end
        %% Use the cooefficients from the training data to predict the sample that was left out
        %     samplename = uniq{k};
        %     goodProbesIdx(:).(samplename) = find(mainA(:,1) ~= 0);
        %
        %     Uhat(testSample(:), :) = X(testSample(:), :) * A;
        %     Vhat(testSample(:), :) = Y(testSample(:), :) * B;
    end
    [modelcorr(:,numpc), pVal(:,numpc)] = diag(corr(Y,yfit));
end
% %% Plot Main PLSR
% 
% mainPslr = figure(1);
% subplot(2,2,1)
% plot(1:length(mainPCTVAR),cumsum(100*mainPCTVAR(2,:)),'-bo');
% xlabel('Number of PLS components');
% ylabel('Percent Variance Explained in y');
% 
% % subplot(2,2,2)
% % yfit = [ones(size(X,1),1) X] * mainBETA;
% % residuals = Y - yfit;
% % stem(residuals)
% % xlabel('Observations');
% % ylabel('Residuals');
% 
% subplot(2,2,3)
% scatter(1:length(vipScore), vipScore,'x')
% hold on
% scatter(indVIPs, vipScore(indVIPs),'rx')
% plot([1 length(vipScore)],[1 1],'--k')
% hold off
% axis tight
% xlabel('Predictor Variables')
% ylabel('Trait Scores')
% 
% subplot(2,2,4)
% plot(0:26,mainMSE(2,:),'b-o');
% xlabel('Number of components');
% ylabel('Estimated Mean Squared Prediction Error');
% legend({'PLSR'},'location','NE');
%% Plot the lasso predictions compared to the actual Comp for a sanity check

MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
    ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];

% check = figure(2);
% subplot(2,2,1)
% modelcorrLOO = corr(Uhat,Vhat);
% heatmap(cornamesV, cornamesU, modelcorrLOO, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
% xlabel('Predicted V (Vhat) PSLR');
% ylabel('Predicted U (Uhat) PSLR');
% 
% subplot(2,2,2)
% modelcorrPSLR = corr(mainU, mainV);
% heatmap(cornamesV, cornamesU, modelcorrPSLR, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
% xlabel('mainV PSLR');
% ylabel('mainU PSLR');

%% %%%%%%%%%% Trait Heatmap Section %%%%%%%%%%%%%%%%

% subplot(2,2,3)
% modelcorrLOO = corr(Y, Vhat);
% heatmap(cornamesV, [Vblnames,cornamesU'], modelcorrLOO, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
% xlabel('Predicted PSLR Components');
% ylabel('Traits');
% 
% subplot(2,2,4)
% modelcorrPSLR = corr(mainY, mainV);
% heatmap(cornamesV, Vblnames, modelcorrPSLR, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
% xlabel('Main PSLR Components');
% ylabel('Traits');

%% predict traits

preditTraits = figure(3);
predictVSactual = corr(mainY, yfit);
heatmap(Vblnames, Vblnames, modelcorr, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Predicted traits');
ylabel('Traits');

for errorall = 1:width(mainY)
    maeall(:,errorall) = mean(abs(mainY(:,errorall) - yfit(:,errorall)));
end
%% Save everything

% % PNGs
% saveas(check, figNameALL);
% saveas(mainPslr, 'PSLRinfo_methyl_PSLRappend.png');
saveas(preditTraits, 'PSLR_predict_vs_actual_traits.png');



% saveas(CC_heatmapLOO, heatmapNameLOO_Traits)
% saveas(CC_heatmapCCA, heatmapNameCCA_Traits)
% saveas(checkV, heatmapNameLOO);
% saveas(checkCCA, heatmapNameCCA);
