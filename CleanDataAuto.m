clear; clc; close all;

%% Go completely automatic?

liltest = input('Use Lillitest Suggestions(Y/N): ', 's');

%% Import traits from metadata

T = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');
T = T(:,3:end);
traits = T.Properties.VariableNames;
traits = string(traits);
Tvalues = table2array(T);

% Preclean the data (remove NaN values)
Tvalues = knnimpute(Tvalues);

%% Pick traits (columns) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
manipulatedVariable = 1:width(T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
valuesDesiredTraits = Tvalues(:, manipulatedVariable);
desiredTraits = traits(manipulatedVariable);
transformed = zeros(height(valuesDesiredTraits), width(manipulatedVariable));

for i = 1:numel(manipulatedVariable)
    %% declare p values and names
p = zeros(10,1);

    names = [" (Raw)", " (^3)", " (^2)", " (^1^/^2)", " (Log(n))", " (Log10)", " (^-^1^/^2)", " (Inverse)", ...
         " (Inverse^2)", " (Inverse^3)"];
    histogramNames = ["(1): Raw", "(2): ^3", "(3): ^2", "(4): ^1^/^2", "(5): Log(n)", ... 
        "(6): Log10", "(7): ^-^1^/^2", "(8): Inverse", "(9): Inverse^2", "(10): Inverse^3"];
     
%% Preclean the data part 2: remove zeros by shifting data by 1% of the range
     
     for m = 1:height(valuesDesiredTraits) 
             valuesDesiredTraits(m, i) = (0.1 * range(valuesDesiredTraits(:,i))) + valuesDesiredTraits(m,i);             
     end
     
    manipulated(:,1) = valuesDesiredTraits(:, i);          % raw
    manipulated(:,2) = valuesDesiredTraits(:, i).^(3);     % cubed
    manipulated(:,3) = valuesDesiredTraits(:, i).^(2);     % squared
    manipulated(:,4) = valuesDesiredTraits(:, i).^(.5);    % sqrt
    manipulated(:,5) = log(valuesDesiredTraits(:, i));     % natural log
    manipulated(:,6) = log10(valuesDesiredTraits(:, i));   % log base 10
    manipulated(:,7) = valuesDesiredTraits(:, i).^(-.5);   % inverse sqrt
    manipulated(:,8) = valuesDesiredTraits(:, i).^(-1);    % inverse
    manipulated(:,9) = valuesDesiredTraits(:, i).^(-2);    % inverse^2
    manipulated(:,10) = valuesDesiredTraits(:, i).^(-3);   % inverse^3
        
    %% lillietest to find 'best' transformation
    for j = 1:width(manipulated)    

       [~,p(j)] = lillietest(manipulated(:,j), 'MCTol', 1e-3); 

    end
    
     for x = 1:height(p)
         if p(x) == 0
             p(x) = inf;
         end
     end
     
    minT = min(p);
    LillieBest = find(abs(p-minT) < 1e-6, 1);
    %% Ask user for the best normal graph
    if liltest == 'N'      
        figure(1)
        
        for k = 1:10
            subplot(3, 4, k);
            histogram(manipulated(:,k));
            xlabel(histogramNames(k));
        end
        
        fprintf('Lillietest thinks graph %d is best\n', LillieBest);
        best = input('Most Normal Graph (1-10): ');
        
        transformed(:,i) = manipulated(:,best); % take the best clean data from lillietest
        desiredTraits(i) = append(desiredTraits(i), names(best)); % rename clean data
    end

end
    
%% Save the best transformed data

cleanData = array2table(transformed, 'VariableNames', desiredTraits);
writetable(cleanData, 'CleanTraitData.csv');
