clear; clc; close;

%% Read in the methylation matrix
methylp = readmatrix('MethMatrixCG-HumanBuccalTurk95_10-100.csv');
% readmatrix DOESNT WORK IF THERE IS INCONSISTENT DATATYPES WITHIN A COLUMN
% Only have doubles if possible, it doesn't read strings
methylp = methylp(:,2:end);

%% Import traits from metadata, and import the sample names we have methylation values for

T = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');

samplenames = sort(textread('Turkish10_100_IDs.txt','%s'));
traitnames = T.No;

%% Compare names from metadata and methylation value names to make sure they are the same

kitid = zeros(length(traitnames),1);

for kit=1:numel(samplenames)
    if length(find(strcmp(traitnames,samplenames{kit}))) == 1
        kitid(kit) = find(strcmp(traitnames,samplenames{kit}));
    end
end
kitid = kitid';

% [v, w] = unique(kitid, 'stable');
% dupidx = setdiff(1:numel(samplenames),w);

%% Filter out names that we don't have methylation or metadata on
idxA = find(kitid > 0);
idxB = kitid(1:95)';
% all valid indexes for finding the right ages
%% Declare X and Y for lasso regression
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

X = methylp(idxA,:);

Y = T.Age(idxB, :);
%Should be all the ages, be sure that the metadata ID's (ProsperIds)
% align with the ID's with methylation data (b.txt)

%% Lasso Regression
i = unique(samplenames);
numUniqueSamples = numel(i);
numsamples = numel(idxA);

for traitidx = 1:1
    traitidx
    for k=1:numUniqueSamples
        k
        idx = [1:numsamples]';
        i = unique(samplenames);
        id = table(idx, samplenames);
        
        j = 1;
        count = 1;
        exit = 0;
        while count <= numsamples
            if isequal(i{k},id.samplenames{j}) == true
                id(j,:) = [];
                if count < numsamples && isequal(i{k},id.samplenames{j})
                    id(j,:) = [];
                end
                exit = 1;
            end
            if exit == 1
                break
            end
            count = count + 1;
            j = j + 1;
        end
        
        XTrain = X(id.idx,:);
        yTrain = Y(id.idx,traitidx);
        
        % L1-constrained linear regression/fits: absolute value (lasso)
        % L2-constrained linear regression/fits: squared (elastic
        % net/ridge)
        [B,FitInfo] = lasso(XTrain,yTrain,'Alpha',1.0,'CV',10);
        % We are doing leave one out cross validation for the samples, i.e.
        % using 99.99% of the data to train the algorythm, then we are
        % testing the algorythm with the one that was left out. We are
        % simultaniously doing 10 fold cross validation (leaving out 10% of
        % the data to test and using 90% of the data to train) for determining
        % the best value of lambda using lasso
        idxLambda1SE = FitInfo.Index1SE;
        % Common convention to always use the lambda value that is 1SE to
        % the left of the minimum lambda value (see MATLAB lasso regression
        % video)
        coef = B(:,idxLambda1SE);
        coef0 = FitInfo.Intercept(idxLambda1SE);
        
        yhat(k,traitidx) = X(k,:)*coef + coef0;
        % ERROR this needs to become an average between the two different
        % methylation rows, this just uses the top row
        
    end
    
end

%% Calculate correlation (variance) between chronological age and predicted age after lasso regression

ageId = unique(idxB);
Y = T.Age(ageId);

modelcorr = corr(yhat,Y);
%corr function determines correlation

%% Plot the lasso predictions compared to the actual age for a sanity check
for z=1:1
    
    check = figure;
    scatter(Y(:,i),yhat(:,i));
    
end
xlabel('Actual Age');
ylabel('Predicted Age');
saveas(check, 'turkLassoTry2.png');
%% Calculate Median Absolute Error

maeage = median(abs(Y(:,1) - yhat(:,1)));

%% Print out the sample names, predicted ages, and actual ages to manually compare

size = [numUniqueSamples, 4];
varTypes = ["string", "double", "double", "double"];
varNames = ["ID", "actualAge", "predictAge", "difference"];
Tpp = table('Size', size, 'VariableTypes', varTypes, 'VariableNames', varNames);
% --> should be all the names and values

Tpp.ID = i(:);
Tpp.actualAge = Y(:,1);
Tpp.predictAge = yhat(:,1);
Tpp.difference = (Tpp.predictAge(:) - Tpp.actualAge(:));

size = [5 2];
varTypes = ["string", "double"];
varNames = ["criticalVariables", "criticalValues"];
critValues = table('Size', size, 'VariableTypes', varTypes, 'VariableNames', varNames);
critValues.criticalVariables = ["MAE"; "coef0"; "Model Correlation"; "idxLambda at 1 SE"; "Sample Size"];
critValues.criticalValues = [maeage; coef0; modelcorr; idxLambda1SE; length(Tpp.actualAge)];

writetable(critValues, 'criticalValuesTurk_10_100.csv');
writetable(Tpp,'lassoTurkresults_10_100.csv');