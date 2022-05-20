%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
 Comparing traits

    Description:

Fit a linear regression model using a matrix input data set.

%}   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close;

%% Read in the phenotypes and the methlyation age predictions
% Use table2array function if dot indexing isnt an option
methylValues = readtable('lassoTurkresults_10_100.csv', 'VariableNamingRule', 'preserve');
IDs = table2array(methylValues(:,1));
methylValues = methylValues(:,2:4);
methylValues = table2array(methylValues);
phenotypesTable = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');
% NOTE: males are labelled as '0' and females are labeled as '1'
phenotypes = phenotypesTable{:,1};
phenotypes = horzcat(phenotypes, phenotypesTable(:,3:end));
phenotypes = table2cell(phenotypes);
actualAge = methylValues(:,1);
predictAge = methylValues(:,2);
phenoData = cell(height(IDs), width(phenotypes));

%% reduce the dataset to only include the metadata for the sample IDs we have

k = 1;
for i=1:height(phenotypes)
    for j = 1:height(IDs)
        if isequal(phenotypes{i,1}, IDs{j})
            phenoData(k,:) = phenotypes(i,:);
            k = k + 1;
        end
    end
end


%% Find the traits with too little information --> this was done manually

 count(1:width(phenoData)) = 0;
 for k = 1:width(phenoData)     
     for i = 1:height(phenoData)
         if isnan(phenoData{i,k})
             count(k) = count(k) + 1;
         end        
     end
 end

%% filter out the markers with not enough information --> this was manually done

 missing = floor(0.25 * height(phenoData));
 
 goodidx = zeros(width(phenoData),1);
 
 for i = 2:width(phenoData)
     if count(i) < missing
         goodidx(i) = i;
     end
 end
 
goodidx = find(goodidx > 0);
%% These are the idices of the columns with sufficient amount of data
 
goodMarkers = cell2mat(phenoData(:, goodidx));
% This is basically phenotypes but only with the markers with sufficient
% data

%% run linear regression on the good data

columnData = zeros(height(actualAge), 1);
multiply = zeros(height(actualAge), 1);
size = [width(goodMarkers), 2];
varTypes = ["cell",  "cell"];
varNames = ["Marker", "Coefficients"];
coeffStore = table('Size', size, 'VariableTypes', varTypes, 'VariableNames', varNames); 
j = 3;

for k = 1:width(goodMarkers)
    coeffStore{k, 1} = phenotypesTable.Properties.VariableNames(j);
    j = j+1;
end

fprintf('Biomarker:\tEstimate\tSE\ttStat\tpValue\n---------------------------------\n');

for k = 1:width(goodMarkers)
    for i = 1:height(actualAge)
        multiply(i,k) = goodMarkers(i, k) .* actualAge(i,1);
        % problem because we aren't using all of the phenotype values -->
        % need to filter out using the indexes from the methyl data
    end
    
    
    X = [goodMarkers(:,k), actualAge, multiply(:,k)];
    mdl = fitlm(X,predictAge);
    
    if mdl.Coefficients.pValue(2) < 0.1
    coeffStore{k, 2} = {mdl.Coefficients};
    
    fprintf('%s: \t%f\t%f\t%f\t%f\n', coeffStore.Marker{k,1}, mdl.Coefficients.pValue(1), ...
        mdl.Coefficients.pValue(2), mdl.Coefficients.pValue(3), mdl.Coefficients.pValue(4));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finding quick correlations


% for k = 1:width(phenotypes)
%     for i = 1:height(actualAge)
%         columnData(i,k) = phenotypes(i, k);
%         multiply(i,k) = phenotypes(i, k) .* actualAge(i,1);
%         % problem because we aren't using all of the phenotype values -->
%         % need to filter out using the indexes from the methyl data
%     end
%     
%     X = [columnData(:,k), actualAge, multiply(:,k)];
%     mdl = fitlm(X,predictAge);
%     
%     if mdl.Coefficients.pValue(2) < 0.1
%     coeffStore{k, 2} = {mdl.Coefficients};
%     end
% end

% hemo_neut_sex = corr(


