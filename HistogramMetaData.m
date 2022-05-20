clear; clc; close;


%% Import traits from metadata
T = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');
traitNames = T.Properties.VariableNames;
traitNames = traitNames(3:end);
Tvalues = table2array(T(:,3:end));

%% Non-linear Transformations
Tvalues(:,4) = (Tvalues(:,4)).^.5; % SpO2
traitNames{4} = "sqrt(SpO2)";

Tvalues(:,5) = (Tvalues(:,5)).^-1; % sytosolic
traitNames{5} = "1/systolic";

Tvalues(:,8) = (Tvalues(:,8)).^2; % neutraphil %
traitNames{8} = "sqrt(neutraphil%)";

Tvalues(:,9) = (Tvalues(:,9)).^.5; % neutraphil #
traitNames{9} = "sqrt(neutraphil#)";

Tvalues(:,10) = (Tvalues(:,10)).^.5; % lymphocyte%
traitNames{10} = "sqrt(lymphocyte%)";

Tvalues(:,12) = (Tvalues(:,12)).^.5; % monocyte%
traitNames{12} = "sqrt(monocyte%)";

Tvalues(:,13) = (Tvalues(:,13)).^.5; % eosinophil%
traitNames{13} = "sqrt(eosinophil%)";

Tvalues(:,14) = (Tvalues(:,14)).^-1; % hemoglobin%
traitNames{14} = "hemoglobin^-^1";

Tvalues(:,15) = (Tvalues(:,15)).^(-3);  % RDW
traitNames{15} = "(RDW)^-^3";

Tvalues(:,17) = log(Tvalues(:, 17)); % nlr
traitNames{17} = "ln(nlr)";

Tvalues(:,18) = (Tvalues(:,18)).^-.5; % Urea
traitNames{18} = "inv.sqrt(Urea)";

Tvalues(:,19) = (Tvalues(:,19)).^-.5; % Creatinine
traitNames{19} = "inv.sqrt(Creatinine)";

Tvalues(:,20) = (Tvalues(:,20)).^-.5; % AST
traitNames{20} = "inv.sqrt(AST)";

Tvalues(:,21) = log10(Tvalues(:, 21)); % ALT
traitNames{21} = "log10(ALT)";

Tvalues(:,22) = log10(Tvalues(:, 22)); % CRP
traitNames{22} = "log10(CRP)";

Tvalues(:,24) = log10(Tvalues(:, 24)); % Cholesterol
traitNames{24} = "log10(Cholesterol)";

Tvalues(:,25) = (Tvalues(:,25)).^-.5; % Triglyceride
traitNames{25} = "inv.sqrt(Triglyceride)";

Tvalues(:,26) = (Tvalues(:,26)).^.5; % LDL
traitNames{26} = "sqrt(LDL)";

Tvalues(:,27) = log10(Tvalues(:, 27)); % VLDL
traitNames{27} = "log10(VLDL)";

Tvalues(:,28) = (Tvalues(:,28)).^-1; % glucose
traitNames{28} = "1/glucose";

Tvalues(:,29) = log10(Tvalues(:, 29)); % Ferritin
traitNames{29} = "log10(Ferritin)";

Tvalues(:,30) = (Tvalues(:,30)).^.5; % CK-MB
traitNames{30} = "sqrt(CK-MB)";

Tvalues(:,31) = log10(Tvalues(:, 31)); % Troponin
traitNames{31} = "log10(Troponin)";

Tvalues(:,32) = log(Tvalues(:, 32)); % D-Dimer
traitNames{32} = "ln(D-Dimer)";

%% Plot histograms

f1 = figure(1);

tagline = 'Histogram_allTraits_Non-linear';
format = '.png';

for i = 1:width(Tvalues)
    subplot(8, 4, i);
    histogram(Tvalues(:,i));
    xlabel(traitNames(i));
    ylabel('frequency');
end

    savename = strcat(tagline, format);
    saveas(f1, savename);