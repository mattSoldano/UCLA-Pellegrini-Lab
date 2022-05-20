clear; clc; close;

%% Read in the methylation matrix
methylT = readtable('MethMatrixCG-HumanBuccalTurk95_10-100.csv', 'VariableNamingRule', 'preserve');
methylCoords = methylT.Properties.VariableNames(2:end);
samplenames = sort(methylT.sampleID);
methyldata = readmatrix('MethMatrixCG-HumanBuccalTurk95_10-100.csv');
methyldata = methyldata(:,2:end);

%% Import traits from metadata, and import the sample names we have methylation values for

T = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');

%% Import traits and import the sample names we have methylation values for

Vblnames = T.Properties.VariableNames;
traitnames = sort(textread('All_Turkish10_100_IDs.txt', '%s'));
% Same IDs as the trait metadata excel file
%% Compare names from metadata and methylation value names to make sure they are the same

kitid = zeros(length(traitnames),1);

for kit=1:numel(samplenames)
    if length(find(strcmp(traitnames,samplenames{kit}))) == 1
        kitid(kit) = find(strcmp(traitnames,samplenames{kit}));
    end
end
kitid = kitid';
methylIdx = find(kitid > 0);
traitIdx = kitid(1:95);

%% Lasso section

X = methyldata(methylIdx,:);

Y = T.Age(traitIdx, :);
% CCs

%% Lasso Regression
i = unique(samplenames);
numUniqueSamples = numel(i);
numsamples = height(Y);
yhat = zeros(numsamples, 1);
IDX = 1;

for k=1:numUniqueSamples
    k
    idx = [1:numsamples]';
    i = unique(samplenames);
    id = table(idx, samplenames);

    j = 1;
    count = 1;
    exit = 0;
    % removing one individual (both methylation and age values) logic for lasso function
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
    yTrain = Y(id.idx, 1);

    if k == 26
        disp('hi');
    end

    [Beta,FitInfo] = lasso(XTrain,yTrain,'Alpha',1.0,'CV',10);
    
    idxLambda1SE = FitInfo.Index1SE;
    
    coef = Beta(:,idxLambda1SE); 
    coef0 = FitInfo.Intercept(idxLambda1SE);
    z = 1;
    hatCount = 1;
    exit = 0;
    
    % yhat predictions for each individual (predicting with both
    % methylation samples)
    while hatCount <= numsamples
        if isequal(i{k},samplenames{z}) == true
            yhat(IDX, 1) = X(IDX,:)*coef + coef0;
            z = z + 1;
            IDX = IDX + 1;
            if hatCount < numsamples && isequal(i{k},samplenames{z})
                yhat(IDX, 1) = X(IDX,:)*coef + coef0;
                IDX = IDX + 1;
            end
            exit = 1;
        end
        if exit == 1
            break
        end
        hatCount = hatCount + 1;
        z = z + 1;
    end
end

%% Calculate correlation (variance) between actual Comp and predicted Comp after lasso regression

modelcorr = corr(yhat,Y);

%% Plot the lasso predictions compared to the actual Comp for a sanity check

for z=1:1
    
    check = figure(3);
    scatter(Y(:),yhat(:));
    
end
xlabel('Actual Age');
ylabel('Predicted Age');
%% Calculate Median Absolute Error

maeComp = median(abs(Y(:,1) - yhat(:,1)));

%% Print out the sample names, predicted Comp, and actual Comp to manually compare

size = [numsamples, 4];
varTypes = ["string", "double", "double", "double"];
varNames = ["ID", "actual", "predict", "difference"];
Tpp = table('Size', size, 'VariableTypes', varTypes, 'VariableNames', varNames);

Tpp.ID = traitnames(traitIdx(:));
Tpp.actual = Y(:,1);
Tpp.predict = yhat(:,1);
Tpp.difference = (Tpp.predict(:) - Tpp.actual(:));

size = [5 2];
varTypes = ["string", "double"];
varNames = ["criticalVariables", "criticalValues"];
critValues = table('Size', size, 'VariableTypes', varTypes, 'VariableNames', varNames);
critValues.criticalVariables = ["MAE"; "coef0"; "Model Correlation"; "idxLambda at 1 SE"; "Sample Size"];
critValues.criticalValues = [maeComp; coef0; modelcorr; idxLambda1SE; length(Tpp.actual)];

%% Save key info
writetable(critValues, 'criticalValuesTurk_AgeTest_raw_age.csv');
writetable(Tpp,'lassoTurkresults_AgeTest_raw_age.csv');
saveas(check, 'turkLasso_AgeTest_raw_age.png'); % scatter plot
