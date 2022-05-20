clear; clc; close;

%% Read in the methylation matrix
methylp = readmatrix('MethMatrixCG-HumanBuccalTurk95_10-100.csv');
% readmatrix DOESNT WORK IF THERE IS INCONSISTENT DATATYPES WITHIN A COLUMN
% Only have doubles if possible, it doesn't read strings
methylp = methylp(:,2:end);

%% Import traits from PC1, and import the sample names we have methylation values for

T = readtable('pcaData_AllPCs.csv', 'VariableNamingRule', 'preserve');
PC6 = T.PC6;

samplenames = sort(textread('Turkish10_100_IDs.txt','%s'));
traitnames = sort(textread('All_Turkish10_100_IDs.txt', '%s'));

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
% all valid indexes for finding the right PC1s
%% Declare X and Y for lasso regression
% Make sure X has the correct amount of samples and that the samples are
% properly aligned

X = methylp(idxA,:);

Y = PC6(idxB, :);
%Should be all the PC1s, be sure that the metadata ID's (ProsperIds)
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

%% Calculate correlation (variance) between chronological PC1 and predicted PC1 after lasso regression

PC6_Id = unique(idxB);
Y = PC6(PC6_Id);

modelcorr = corr(yhat,Y);
%corr function determines correlation

%% Plot the lasso predictions compared to the actual PC1 for a sanity check

for z=1:1
    
    check = figure(1);
    scatter(Y(:),yhat(:));
    
end
xlabel('Actual PC6');
ylabel('Predicted PC6');
%% Calculate Median Absolute Error

maePC6 = median(abs(Y(:,1) - yhat(:,1)));

%% Print out the sample names, predicted PC1s, and actual PC1s to manually compare

size = [numUniqueSamples, 4];
varTypes = ["string", "double", "double", "double"];
varNames = ["ID", "actualPC6", "predictPC6", "difference"];
Tpp = table('Size', size, 'VariableTypes', varTypes, 'VariableNames', varNames);
% --> should be all the names and values

Tpp.ID = i(:);
Tpp.actualPC6 = Y(:,1);
Tpp.predictPC6 = yhat(:,1);
Tpp.difference = (Tpp.predictPC6(:) - Tpp.actualPC6(:));

size = [5 2];
varTypes = ["string", "double"];
varNames = ["criticalVariables", "criticalValues"];
critValues = table('Size', size, 'VariableTypes', varTypes, 'VariableNames', varNames);
critValues.criticalVariables = ["MAE"; "coef0"; "Model Correlation"; "idxLambda at 1 SE"; "Sample Size"];
critValues.criticalValues = [maePC6; coef0; modelcorr; idxLambda1SE; length(Tpp.actualPC6)];

%% Save key info
writetable(critValues, 'criticalValuesTurk_PC6.csv');
writetable(Tpp,'lassoTurkresults_PC6.csv');
saveas(check, 'turkLassoPC6.png'); % scatter plot