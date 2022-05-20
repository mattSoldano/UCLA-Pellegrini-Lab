clear; clc;

% DO NOT INCLUDE IN GITHUB

%% Make a string with chromosomes and positions seperated
data = readcell('PC_Trait_AssociationModel_SiteSpecific.xlsx', 'VariableNamingRule', 'preserve');
Vblnames = data(1,1:end-1);
data = data(2:end,1:end-1);
foundSites = data(:,1);

for i=1:numel(foundSites)
    foundSites{i} = convertCharsToStrings(foundSites{i});
end

for i=1:numel(foundSites)
    CHRs{i,1} = extractBefore(foundSites{i}, ':');
    positions{i,1} = str2double(extractAfter(foundSites{i}, ':'));
end

Sites = [CHRs, positions];

%% Import HG38 database to find gene associations

database = readcell('SiteAssociationGeneCheck.csv');
GeneVblnames = database(1,1:end-2);
database = database(2:end,1:end-2);
checkStartEnd = database(:,3:4);
checkChr = database(:,2);

for i=1:numel(checkChr)
    checkChr{i} = convertCharsToStrings(checkChr{i});
end

geneCheck = [checkChr, checkStartEnd]; 

%% Find regions of genes and their associations with CpGs

for s = 1:height(Sites)
   
    if strcmp(Sites{s,1}, geneCheck{s,1}) 
        % see video from 5/16 meeting with Pellegrini
        
    end
    
end

% sites = readcell('human_probes.txt');
% sites = sites(:,2:end);
% for i=1:height(sites)
%     if isnumeric(sites{i,1})
%         NumSitesIdx(i,1) = i;
%     end
% end
% 
% NumSitesIdx = find(NumSitesIdx ~= 0);
% 
% for i=1:numel(NumSitesIdx)
%     sites{i,1} = num2str(sites{i,1});
% end
% for i=1:height(sites)
%     sites{i,1} = strcat('chr', sites{i,1});
% end
% writecell(sites, 'probeSites.txt');