
pause
clear; clc; close;

%{ 
Try to just use beta values from turk data on proteome data to see if the 
beta values are robust
%}
%% Read in the methylation matrices --> rename sampleIDs (manually in a txt file unsorted)
cd ~/Pellegrini/firstTurkDataset/Proteomics_Data
methylT = readtable('MethMatrixCG-HumanBloodFitness96manual_20-100.csv', 'VariableNamingRule', 'preserve');
methylData = table2array(methylT(:, 2:end));
methylSites = methylT.Properties.VariableNames(2:end);
samplenames = methylT.sampleId;

%% make a single methylation matrix using all the samples
%% Import traits and import the sample names we have methylation values for

T = readtable('Traits_ProsperProteome103_120821mjt.csv', 'VariableNamingRule', 'preserve');
Traits = T.Properties.VariableNames(2:end);
traitdata = table2array(T(:, 2:end));
traitnames = T.sampleID;

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

load('TurkBeta');
% load('Vblnames');
% load('allnames');
sameTraitsTurk = [1, 2, 3, 10, 16, 17, 18, 19, 22, 23];
sameTraitsProteome = [3, 2, 53, 54, 79, 80, 91, 83, 85, 84];

% 56 and 57 neuotrophils and lymphocytes
Traits = Traits(sameTraitsProteome);

X = methylData;
Y = traitdata(traitIdx, :);
Y = Y(:, sameTraitsProteome);
BETA = TurkBeta(:, sameTraitsTurk);
% were different probesets used? Why are there a different number of probes?
yfit = [ones(numsamples,1) X] * BETA;

MYmap = [linspace(0,1,25)', linspace(0,1,25)', linspace(1,1,25)' ...
    ; linspace(1,1,25)' linspace(1,0,25)' linspace(1,0,25)'];

%% predict traits

for t = 1:width(mainU)
    cornamesU(t,1) = sprintf("U%d", t);
    cornamesV(t,1) = sprintf("V%d", t);
    cornames(t,1) = sprintf("Comp. %d", t);
end

preditTraits1 = figure(1);
predictVSactual = corr(yfit, Y);
heatmap(vbls_comps, vbls_comps, predictVSactual, 'ColorLimits', [-1, 1], 'Colormap', MYmap);
xlabel('Predicted traits');
ylabel('Traits');

for errorall = 1:width(Ycat)
    maeall(:,errorall) = mean(abs(Ycat(:,errorall) - yfit(:,errorall)));
end