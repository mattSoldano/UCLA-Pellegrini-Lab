%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
 PCA traits for Metabolic Syndrom

    Description:
%}   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close;

%% Step 1: data standardizion
T = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');
data = readmatrix('Turkish_Buccal_Traits.csv');
dataCellType = data(:,25:30);

ID = readcell('Turkish_Buccal_Traits.csv');
ID = ID(2:end,1);

variableAvgs = mean(dataCellType);
variableSDs = std(dataCellType);
nData = normalize(dataCellType, variableAvgs, variableSDs);

% How should I fill in the data if there is NaN values? Normalizing data
% won't work with NaN values

% Remove NaN values for normData

for i = 1:width(nData)
    for k = 1:height(nData)
        if isnan(nData(k,i)) == false
        normData(k,i) = nData(k,i);
        end
    end
end

%% Step 2: covariance matrix computation
covM = (normData' * normData) / height(normData);

% Remove NaN values for covariance

for i = 1:width(covM)
    for k = 1:height(covM)
        if isnan(covM(k,i)) == false
        covarianceM(k,i) = covM(k,i);
        end
    end
end

%% Step 3: finding and sorting the eigenvalues and eigenvectors based on the 
% absolute value of the eigenvalues
[coeffOrth, eValues] = eig(covarianceM, 'vector');
[~, idx] = sort(abs(eValues), 'descend');
coeffOrth = coeffOrth(:,idx);

%% Step 4: projecting the data into the new space whose bases are the 
% eigenvectors corresponding to the top-k values
pcaData = normData * coeffOrth;

%% Step 5: Run biPlot

% save variables from covidCountries table
vbls = T.Properties.VariableNames;
Tvbls = vbls(:, 25:30);
% take only the first two components
E2D = coeffOrth(:,1:2);
pcaData2D = pcaData(:,1:2);
% Run biplot
PCAplot = figure;
biplot(E2D, 'Scores', pcaData2D, 'Varlabels', Tvbls);
saveas(PCAplot, 'PCA_biplot_Metabolic.png');


function normalizedData = normalize(dataset, avgs, SDs)
normalizedData = zeros(height(dataset), width(dataset));
for i = 1:width(dataset)
    for j = 1:height(dataset)
        normalizedData(j,i) = ((dataset(j,i) - avgs(i)) ./ SDs(i));
    end
end
end