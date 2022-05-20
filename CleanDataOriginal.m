clear; clc; close;

%% Import traits from metadata

T = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');
T = T(:,3:end);
traits = T.Properties.VariableNames;
Tvalues = table2array(T);
Tvalues = knnimpute(Tvalues);

%% Pick trait

manipulatedVariable = 12;

%% change values of zero to small values (max 10% of min) 
     n = 1;
     for m = 1:height(Tvalues)
        if Tvalues(m, manipulatedVariable) ~= 0
             minV(n) = Tvalues(m, manipulatedVariable);
             n = n + 1;
        end
     end
     
     minV = min(minV);
     
     for z = 1:height(Tvalues)
         if Tvalues(z, manipulatedVariable) == 0
             Tvalues(z, manipulatedVariable) = 0.1 * minV;
         end
     end

%% Make the manipulated array
manipulated = zeros(height(Tvalues), 10);
names = ["raw", "cube", "square", "sqrt", "nlog", "log10", "-sqrt", "inverse", ...
    "inverse^2", "inverse^3"];

    manipulated(:,1) = Tvalues(:, manipulatedVariable);
    manipulated(:,2) = Tvalues(:, manipulatedVariable).^(3);
    manipulated(:,3) = Tvalues(:, manipulatedVariable).^(2);
    manipulated(:,4) = Tvalues(:, manipulatedVariable).^(.5);
    manipulated(:,5) = log(Tvalues(:, manipulatedVariable));
    manipulated(:,6) = log10(Tvalues(:, manipulatedVariable));
    manipulated(:,7) = Tvalues(:, manipulatedVariable).^(-.5);
    manipulated(:,8) = Tvalues(:, manipulatedVariable).^(-1);
    manipulated(:,9) = Tvalues(:, manipulatedVariable).^(-2);
    manipulated(:,10) = Tvalues(:, manipulatedVariable).^(-3);

    figname = traits{manipulatedVariable};
    
figure(1) 

disp(figname)
disp(manipulatedVariable)

for k = 1:10
subplot(3, 4, k);
histogram(manipulated(:,k));
xlabel(names(k));
end
