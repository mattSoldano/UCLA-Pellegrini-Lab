clear; clc; close;

%% Step 1: data standardizion
T = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');
data = readmatrix('Turkish_Buccal_Traits.csv');
dataCellType = data(:,9:19);
dataCellType = horzcat(dataCellType,data(:,24));
% save variables
vbls = T.Properties.VariableNames;
Tvbls = vbls(:, 9:19);
Tvbls = horzcat(Tvbls,vbls(:,24));

%% Non-linear Transformations
% dataCellType(:,4) = (dataCellType(:,4)).^.5; % SpO2
% Tvbls{4} = 'sqrt(SpO2)';
% 
% dataCellType(:,5) = (dataCellType(:,5)).^-1; % sytosolic
% Tvbls{5} = '1/systolic';
% 
dataCellType(:,2) = (dataCellType(:,2)).^2; % neutraphil %
Tvbls{2} = 'sqrt(neutraphil%)';

dataCellType(:,3) = (dataCellType(:,3)).^.5; % neutraphil #
Tvbls{3} = 'sqrt(neutraphil#)';

dataCellType(:,4) = (dataCellType(:,4)).^.5; % lymphocyte%
Tvbls{4} = 'sqrt(lymphocyte%)';

dataCellType(:,6) = (dataCellType(:,6)).^.5; % monocyte%
Tvbls{6} = 'sqrt(monocyte%)';

dataCellType(:,7) = (dataCellType(:,7)).^.5; % eosinophil%
Tvbls{7} = 'sqrt(eosinophil%)';

dataCellType(:,8) = (dataCellType(:,8)).^-1; % hemoglobin%
Tvbls{8} = 'hemoglobin^-^1';

dataCellType(:,9) = (dataCellType(:,9)).^(-3);  % RDW
Tvbls{9} = '(RDW)^-^3';

dataCellType(:,11) = log(dataCellType(:, 11)); % nlr
Tvbls{11} = 'ln(nlr)';

dataCellType(:,12) = log10(dataCellType(:, 12)); % CRP
Tvbls{12} = 'log10(CRP)';

% dataCellType(:,24) = log10(dataCellType(:, 24)); % Cholesterol
% Tvbls{24} = "log10(Cholesterol)";

% dataCellType(:,25) = (dataCellType(:,25)).^-.5; % Triglyceride
% Tvbls{25} = "inv.sqrt(Triglyceride)";
% 
% dataCellType(:,26) = (dataCellType(:,26)).^.5; % LDL
% Tvbls{26} = "sqrt(LDL)";
% 
% dataCellType(:,27) = log10(dataCellType(:, 27)); % VLDL
% Tvbls{27} = "log10(VLDL)";
% 
% dataCellType(:,28) = (dataCellType(:,28)).^-1; % glucose
% Tvbls{28} = "1/glucose";
% 
% dataCellType(:,29) = log10(dataCellType(:, 29)); % Ferritin
% Tvbls{29} = "log10(Ferritin)";
% 
% dataCellType(:,30) = (dataCellType(:,30)).^.5; % CK-MB
% Tvbls{30} = "sqrt(CK-MB)";
% 
% dataCellType(:,31) = log10(dataCellType(:, 31)); % Troponin
% Tvbls{31} = "log10(Troponin)";
% 
% dataCellType(:,32) = log(dataCellType(:, 32)); % D-Dimer
% Tvbls{32} = "ln(D-Dimer)";

%% histogram

f1 = figure(1);
title('Immune Cell Types Variance')
tagline = 'Histogram';
seperator = '_';
format = '.png';

for i = 1:width(dataCellType)
    subplot(4, 3, i);
    histogram(dataCellType(:,i));
    xlabel(Tvbls(i));
    ylabel('Occurances');
end

    savename = strcat(tagline, '_immune_cell_types', format);
    saveas(f1, savename);