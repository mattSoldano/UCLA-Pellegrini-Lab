clear; clc; close;

allTraits = readtable('Turkish_Buccal_Traits.csv', 'VariableNamingRule', 'preserve');
traitNames = allTraits.Properties.VariableNames;
traitNames = horzcat(traitNames(9:19), traitNames(24));
traits = horzcat(table2array(allTraits(:,9:19)), ...
    table2array(allTraits(:,24)));


f1 = figure(1);
title('Immune Cell Types Variance')
tagline = 'Histogram';
seperator = '_';
format = '.png';

for i = 1:width(traits)
    subplot(4, 3, i);
    histogram(traits(:,i));
    xlabel(traitNames(i));
    ylabel('Occurances');

end

savename = strcat(tagline, '_immune_cell_types', format);
    saveas(f1, savename);