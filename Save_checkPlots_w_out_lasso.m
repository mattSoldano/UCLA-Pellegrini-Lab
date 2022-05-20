clear; clc; close;

T = readtable('lassoTurkResults_PC4_Transformed.csv');

Y = T.actualPC4;
yhat = T.predictPC4;

for z=1:1
    
    check = figure(1);
    scatter(Y(:),yhat(:));
    
end
xlabel('Actual PC4');
ylabel('Predicted PC4');
saveas(check, 'turkLassoPC4_Transformed.png'); % scatter plot