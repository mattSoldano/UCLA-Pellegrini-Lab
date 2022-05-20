%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
 PCA traits for Cell Type

    Description:
%}   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close;

%% Step 1: data standardizion
T = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');
data = readmatrix('Turkish_Buccal_Traits.csv');
dataCellType = data(:,9:19);
dataCellType = horzcat(dataCellType,data(:,24));

ID = readcell('Turkish_Buccal_Traits.csv');
ID = ID(2:end,1);

variableAvgs = mean(dataCellType);
variableSDs = std(dataCellType);
normData = normalize(dataCellType, variableAvgs, variableSDs);

% % How should I fill in the data if there is NaN values? Normalizing data
% % won't work with NaN values, so I just removed them

% Remove NaN values for normData --> would be what to do if there is NaN
% values for the data we are interested in (don't have a script to fix the
% variable names when NaNs are removed)

% for i = 1:width(nData)
%     for k = 1:height(nData)
%         if isnan(nData(k,i)) == false
%         normData(k,i) = nData(k,i);
%         end
%     end
% end

%% Step 2: covariance matrix computation
covarianceM = (normData' * normData) / height(normData);

% % Remove NaN values for covariance

% for i = 1:width(covM)
%     for k = 1:height(covM)
%         if isnan(covM(k,i)) == false
%         covarianceM(k,i) = covM(k,i);
%         end
%     end
% end

%% Step 3: finding and sorting the eigenvalues and eigenvectors based on the 
% absolute value of the eigenvalues
[coeffOrth, eValues] = eig(covarianceM, 'vector');
[~, idx] = sort(abs(eValues), 'descend');
coeffOrth = coeffOrth(:,idx);

%% Step 4: projecting the data into the new space whose bases are the 
% eigenvectors corresponding to the top-k values
pcaData = normData * coeffOrth;

%% Step 5: Run biPlot

% save variables
vbls = T.Properties.VariableNames;
Tvbls = vbls(:, 9:19);
Tvbls = horzcat(Tvbls,vbls(:,24));

% Pick the components to observe:

E2D = coeffOrth(:,5:6);
pcaData2D = pcaData(:,5:6);

% Run biplot
PCAplot = figure;
biplot(E2D, 'Scores', pcaData2D, 'Varlabels', Tvbls);

% Save Name for Plot:
saveas(PCAplot, 'PCA_biplot_cellType_PC5_6.png');

variableNames = {'PC3', 'PC4'};
pcaData2D = array2table(pcaData2D, 'VariableNames', variableNames);

% Save Name for Data Table:
% writetable(pcaData2D, 'pcaDataPlot_PC1_PC2.csv');

function normalizedData = normalize(dataset, avgs, SDs)
normalizedData = zeros(height(dataset), width(dataset));
for i = 1:width(dataset)
    for j = 1:height(dataset)
        normalizedData(j,i) = ((dataset(j,i) - avgs(i)) ./ SDs(i));
    end
end
end