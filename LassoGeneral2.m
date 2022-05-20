clear; clc; close;

%% Read in the methylation matrix
methylT = readtable('MethMatrixCG-HumanBuccalTurk95_10-100.csv', 'VariableNamingRule', 'preserve');
methylCoords = methylT.Properties.VariableNames(2:end);
samplenames = sort(methylT.sampleID);
methyldata = readmatrix('MethMatrixCG-HumanBuccalTurk95_10-100.csv');
methyldata = methyldata(:,2:end);

%% Import Cannonical Correlation or whatever predictor traits you want

T = readtable('CleanTraitData_ccorrComponents_Rescaled.csv', 'VariableNamingRule', 'preserve');
Vblnames = T.Properties.VariableNames;
% Same IDs as the trait metadata excel file
data = readmatrix('CleanTraitData_ccorrComponents_Rescaled.csv');

%% Pick the component you want to be the predictor

Vblnames = Vblnames(7);
data = data(:,7);

X = methyldata;

Y = data;
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
%corr function determines correlation

%% Plot the lasso predictions compared to the actual Comp for a sanity check

check = figure(3);
scatter(Y(:),yhat(:));

xlabel('Actual Comp7');
ylabel('Predicted Comp7');

% %% Fitline:
% coeffs = polyfit(Y, yhat, 1);
% 
% % Get fitted values
% fittedY = linspace(min(Y), max(Y), 200);
% fittedYhat = polyval(coeffs, fittedY);
% hold on
% 
% plot(fittedY, fittedYhat, 'r-', 'LineWidth', 2);
%% Calculate Median Absolute Error

maeComp = median(abs(Y(:,1) - yhat(:,1)));

%% Print out the sample names, predicted Comp, and actual Comp to manually compare

size = [numsamples, 4];
varTypes = ["string", "double", "double", "double"];
varNames = ["ID", "actual", "predict", "difference"];
Tpp = table('Size', size, 'VariableTypes', varTypes, 'VariableNames', varNames);
% --> should be all the names and values

Tpp.ID = samplenames(:);
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
writetable(critValues, 'criticalValuesTurk_Comp7_NewLogic.csv');
writetable(Tpp,'lassoTurkresults_Comp7_NewLogic.csv');
saveas(check, 'turkLassoComp7_NewLogic.png'); % scatter plot
