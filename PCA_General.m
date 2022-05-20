%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Runs PCA and exports values for the PCs

    Description:
%}   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all;

%% Step 1: data standardizion

% Pick your data table:
metaData = "Turkish_Buccal_Traits.csv";
% Pick the columns from this table:
columns = 2 + [1 2 7:28 30 31];
% Name the save file:
exportName = "pcaData_Raw.csv";
pngName = "PCA_biplot_Raw_PC";

T = readtable(metaData, 'VariableNamingRule', 'preserve');
data = readmatrix(metaData);
data = knnimpute(data(:,columns));
% save variables
Tvbls = T.Properties.VariableNames(columns);

%% Get IDs and normalize data

ID = readcell(metaData);
variableAvgs = mean(data);
variableSDs = std(data);
normData = normalize(data, variableAvgs, variableSDs);

%% Step 2: covariance matrix computation
covarianceM = (normData' * normData) / height(normData);

%% Step 3: finding and sorting the eigenvalues and eigenvectors based on the 
% absolute value of the eigenvalues
[coeffOrth, eValues] = eig(covarianceM, 'vector');
[~, idx] = sort(abs(eValues), 'descend');
coeffOrth = coeffOrth(:,idx);

%% Step 4: projecting the data into the new space whose bases are the 
% eigenvectors corresponding to the top-k values
pcaData = normData * coeffOrth;

%% Step 5: Run biPlot for all PCs
j = 1;
k = 2;
variableNames = cell(1, width(pcaData));

for i = 1:width(pcaData)/2
    E2D = coeffOrth(:,j:k);
    pcaData2D = pcaData(:,j:k);
    
    j = num2str(j);
    k = num2str(k);
    
    labelx = append('Component ', j);
    labely = append('Component ', k);
    
    % Run biplot
    PCAplot = figure;
    biplot(E2D, 'Scores', pcaData2D, 'Varlabels', Tvbls);
    xlabel(labelx);
    ylabel(labely);
    
    % Save Name for Plot:
    
    pngSave = append(pngName, j, '_', k, '.png');
    saveas(PCAplot, pngSave);
    pc = {'PC'};
    jPC = strcat(pc, j);
    kPC = strcat(pc, k);    
    j = str2double(j);
    k = str2double(k);
    variableNames(j) = jPC;
    variableNames(k) = kPC;
    j = j + 2;
    k = k + 2;
    
end

% variableNames = {'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', ...
%     'PC9', 'PC10', 'PC11', 'PC12'};
pcaData = array2table(pcaData, 'VariableNames', variableNames);

writetable(pcaData, exportName);


function normalizedData = normalize(dataset, avgs, SDs)
normalizedData = zeros(height(dataset), width(dataset));
for i = 1:width(dataset)
    for j = 1:height(dataset)
        normalizedData(j,i) = ((dataset(j,i) - avgs(i)) ./ SDs(i));
    end
end
end


%% Non-linear Transformations, no longer needed
% dataCellType(:,4) = (dataCellType(:,4)).^.5; % SpO2
% Tvbls{4} = 'sqrt(SpO2)';
% 
% dataCellType(:,5) = (dataCellType(:,5)).^-1; % sytosolic
% Tvbls{5} = '1/systolic';
% 
% dataCellType(:,2) = (dataCellType(:,2)).^2; % neutraphil %
% Tvbls{2} = 'sqrt(neutraphil%)';
% 
% dataCellType(:,3) = (dataCellType(:,3)).^.5; % neutraphil #
% Tvbls{3} = 'sqrt(neutraphil#)';
% 
% dataCellType(:,4) = (dataCellType(:,4)).^.5; % lymphocyte%
% Tvbls{4} = 'sqrt(lymphocyte%)';
% 
% dataCellType(:,6) = (dataCellType(:,6)).^.5; % monocyte%
% Tvbls{6} = 'sqrt(monocyte%)';
% 
% dataCellType(:,7) = (dataCellType(:,7)).^.5; % eosinophil%
% Tvbls{7} = 'sqrt(eosinophil%)';
% 
% dataCellType(:,8) = (dataCellType(:,8)).^-1; % hemoglobin%
% Tvbls{8} = 'hemoglobin^-^1';
% 
% dataCellType(:,9) = (dataCellType(:,9)).^(-3);  % RDW
% Tvbls{9} = '(RDW)^-^3';
% 
% dataCellType(:,11) = log(dataCellType(:, 11)); % nlr
% Tvbls{11} = 'ln(nlr)';
% 
% dataCellType(:,12) = log10(dataCellType(:, 12)); % CRP
% Tvbls{12} = 'log10(CRP)';
% 
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