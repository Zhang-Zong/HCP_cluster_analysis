function hcp_clus_initialize
% loads HCP data for the clsuter study. Data pulled from the S900 release,
% (sorry exact date not available, was pulled by the
% informatics group). 
% This function looks in the specified folder (basedir) and loops through
% for the tstat files, in the format I saved them in. For each modality, we
% have a cope and tstat from a selected contrast. 
%
% Modality:    Contrast num   Contrast descriptin
% 1. EMOTION       3           FACES-SHAPES
% 2. GAMBLING      6           REWARD-PUNISH
% 3. LANGUAGE      4           STORY-MATH
% 4. RELATIONAL    4           REL-MATCH
% 5. SOCIAL        6           TOM-RANDOM
% 6. WM            11          2BK-0BK
%
% Note that we pulled the analysis combining the left gradient and right
% gradients
% We loop through subject folders, load the tstat file into a matrix, and
% exclude any participants who do not have data for all 6 contrasts of
% interest. 
% Each subject has a folder for each of the above 6 modalities. 
%
% The final data variable is size 822 for the 822 partiicpants with all six
% modalities, 58997 for the verticied on the cortical surface with data,
% and 6 for modalities, ordered as the table above. %
% So data (:,:,4) is a subjects x verticies matrix of tstats for the RELATIONAL data

curdir=pwd; 

% The modalities
mods = {
    'EMOTION'
    'GAMBLING'
    'LANGUAGE'
    'RELATIONAL'
    'SOCIAL'
    'WM'}

%Initalization variables. 
basedir='F:\HCP900/data8/'; %this is where I stashed my copies of the t-maps
cd(basedir)
ids=dir('*')
tfile = 'tstat1.dtseries.nii'
inc(1:899) = 1; %there is 899 subjects in folder

% initializing this saves a lot of trouble, and some computational time
tdat = zeros(899, 96854, 6); 

% this loop goes through the HCP folder for each participant and extracts
% data from all six tasks. Only if all six tasks are available does it
% retain a participant. Partiicpants without all tasks are set inc=0
% start at 3 becuase ids 1 and 2 are not folders ("." and "..")
for pdx = 3:length(ids)
    n=pdx-2 %participant number in matrix)
    ID=ids(pdx).name; %folder)
    name{n}=ID; %generate a list of sub IDs, as cells
   
    %Loop through modalities
    for mdx = 1:6
        modality = mods{mdx};
        
        cd([basedir ID '/' modality ]);
        %if there is a tstat file, read data. Otherwise add zeros
        if exist(tfile, 'file')
             cift = ft_read_cifti(tfile);
             tdat(n,:, mdx) = cift.dtseries;
         else
             tdat(n,:,mdx) = 0;
             inc(n) = 0;
        end
    end
end

%cortex only, we are ignoring the sub cortex for now, it seems to add a lot
%of variance to the data
data=tdat(:,1:64569,:);

% exclude participants who are missing a scan in any modality, and remove
% NaN voxels
data = data(inc==1,~isnan(data(1,:,1)),:);

% lsit of names included in final analysis
name_inc = name(inc==1);


% save base data at this stage. 
%go back to the durectory where it all began. 
cd(curdir)
 % clear useless variables. 
clear ID ans  basedir cift   modality pdx  mdx n ids tdat  tfile
save -v7.3 hcp_data8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize hierarchical clustering
% using full sample

% Basic clustering variables and parameters

% linkage functions
Ze= linkage(data(:,:,1), 'ward');
Zg= linkage(data(:,:,2), 'ward');
Zl= linkage(data(:,:,3), 'ward');
Zr= linkage(data(:,:,4), 'ward');
Zs= linkage(data(:,:,5), 'ward');
Zw= linkage(data(:,:,6), 'ward');

% plot dengrograms to generat "p", which is order of dendrograms
figure; [h t pe] = dendrogram(Ze, 0,'ColorThreshold', 22000);
figure; [h t pg] = dendrogram(Zg, 0,'ColorThreshold', 22000);
figure; [h t pl] = dendrogram(Zl, 0,'ColorThreshold', 22000);
figure; [h t pr] = dendrogram(Zr, 0,'ColorThreshold', 22000);
figure; [h t ps] = dendrogram(Zs, 0,'ColorThreshold', 22000);
figure; [h t pw] = dendrogram(Zw, 0,'ColorThreshold', 22000);
close all % close all those messy graphs

% make a loopaple variable of linkages
Zall(:,:,1)=Ze; 
Zall(:,:,2)=Zg;
Zall(:,:,3)=Zl;
Zall(:,:,4)=Zr;
Zall(:,:,5)=Zs;
Zall(:,:,6)=Zw;

% make a loopable variable of dendrgram order
pall = [ pe;   pg ;  pl;   pr;   ps;   pw ];

% generate distance matricies
for mdx = 1:6
    pdis = pdist(data(:,:,mdx), 'euclidean');
    sdis(:,:,mdx) = squareform(pdis); 
end

% note the data variable will not be saved here
save hcp_hiearechical
