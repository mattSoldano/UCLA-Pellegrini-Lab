clear; clc; close;

%% Import traits from metadata

T = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');
traits = T.Properties.VariableNames(3:end);
Tvalues = table2array(T(:,3:end));
Tvalues = knnimpute(Tvalues);

%% Non-linear Transformations
Tvalues(:,4) = (Tvalues(:,4)).^.5; % SpO2
traits{4} = "sqrt(SpO2)";

Tvalues(:,5) = (Tvalues(:,5)).^-1; % sytosolic
traits{5} = "1/systolic";

Tvalues(:,8) = (Tvalues(:,8)).^2; % neutraphil %
traits{8} = "sqrt(neutraphil%)";

Tvalues(:,9) = (Tvalues(:,9)).^.5; % neutraphil #
traits{9} = "sqrt(neutraphil#)";

Tvalues(:,10) = (Tvalues(:,10)).^.5; % lymphocyte%
traits{10} = "sqrt(lymphocyte%)";

Tvalues(:,12) = (Tvalues(:,12)).^.5; % monocyte%
traits{12} = "sqrt(monocyte%)";

Tvalues(:,13) = (Tvalues(:,13)).^.5; % eosinophil%
traits{13} = "sqrt(eosinophil%)";

Tvalues(:,14) = (Tvalues(:,14)).^-1; % hemoglobin%
traits{14} = "hemoglobin^-^1";

Tvalues(:,15) = (Tvalues(:,15)).^(-3);  % RDW
traits{15} = "(RDW)^-^3";

Tvalues(:,17) = log(Tvalues(:, 17)); % nlr
traits{17} = "ln(nlr)";

Tvalues(:,18) = (Tvalues(:,18)).^-.5; % Urea
traits{18} = "inv.sqrt(Urea)";

Tvalues(:,19) = (Tvalues(:,19)).^-.5; % Creatinine
traits{19} = "inv.sqrt(Creatinine)";

Tvalues(:,20) = (Tvalues(:,20)).^-.5; % AST
traits{20} = "inv.sqrt(AST)";

Tvalues(:,21) = log10(Tvalues(:, 21)); % ALT
traits{21} = "log10(ALT)";

Tvalues(:,22) = log10(Tvalues(:, 22)); % CRP
traits{22} = "log10(CRP)";

Tvalues(:,24) = log10(Tvalues(:, 24)); % Cholesterol
traits{24} = "log10(Cholesterol)";

Tvalues(:,25) = (Tvalues(:,25)).^-.5; % Triglyceride
traits{25} = "inv.sqrt(Triglyceride)";

Tvalues(:,26) = (Tvalues(:,26)).^.5; % LDL
traits{26} = "sqrt(LDL)";

Tvalues(:,27) = log10(Tvalues(:, 27)); % VLDL
traits{27} = "log10(VLDL)";

Tvalues(:,28) = (Tvalues(:,28)).^-1; % glucose
traits{28} = "1/glucose";

Tvalues(:,29) = log10(Tvalues(:, 29)); % Ferritin
traits{29} = "log10(Ferritin)";

Tvalues(:,30) = (Tvalues(:,30)).^.5; % CK-MB
traits{30} = "sqrt(CK-MB)";

Tvalues(:,31) = log10(Tvalues(:, 31)); % Troponin
traits{31} = "log10(Troponin)";

Tvalues(:,32) = log(Tvalues(:, 32)); % D-Dimer
traits{32} = "ln(D-Dimer)";

%% Import PCs
PCtable = readtable('pcaDataTransformed.csv', 'VariableNamingRule', 'preserve');
PCnames = PCtable.Properties.VariableNames;
PCs = table2array(PCtable);

%% Declare X and Y for lasso regression
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

modelcorr = zeros(width(traits), width(PCnames));

for pc = 1:width(PCs)
    for t = 1:width(Tvalues)
        
        Y = Tvalues(:,t);
        X = PCs(:, pc);
        
        %% Determining Correlation

        modelcorr(t, pc) = corr(X,Y);

    end
end

myheatmap = figure(1);

heatmapPCs = heatmap(PCnames, traits, modelcorr);
xlabel('Principal Components');
ylabel('Traits');

saveas(myheatmap, 'PC_Heatmap_NonLinearCorrection.png')