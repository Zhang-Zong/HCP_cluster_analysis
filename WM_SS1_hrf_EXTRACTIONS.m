function WM_SS1_hrf_EXTRACTIONS
% this function examined the hemodynamic response to the 0Bk and 2Bk
% conditions in the working memroy task. 
% two sets of 100 participants were selected based on their score on the
% first principal component of the boostrap agreement matrix for the k=4 WM
% task. These are partiicpants who fall on the uppper or lower extremes of
% those axes, and are therefore the 'strongest activators' or 'most
% deactivators' according to that axis (see regression analyssi of
% component 1, which represents a sort of 'global braina ctivity'
% component). I think
%
% Colin Hawco, Jan 2019. This is hack code, not elegant. 

basedir='F:\HCP900\WM_SS1_atlas_dtseries\'

cd(basedir)

% load the names lsit. n1 is a list of subject ids who have low values
% on the component, n2 are those with high values.
load WM_ss1_ends_namelist.mat

% block duration
dur = floor(27.5/0.72)+4;

% I downlaoded the WM dtseries files (files labeled atlas) from the HCP directories
% this data is unsmoothed. In oder to match the data used in the GLMs I
% will need to smooth 8mm. Also, for convience, I will be concatinating the
% files (LR and RL are separated).

for pdx = 1:length(n1)
    try
        cd(basedir)
        dat(:,:,pdx)=csvread([n1{pdx} '_WM_merged_atlas_glasserROIs.csv']);
        
        % get even onset data
        cd([basedir n1{pdx} '\EVs_LR'])
        
        ev_0bk=[]; ev_2bk=[];
        t=textread('0bk_body.txt');
        ev_0bk=[ev_0bk t(1)];
        t=textread('0bk_tools.txt ');
        ev_0bk=[ev_0bk t(1)];
        t=textread('0bk_faces.txt');
        ev_0bk=[ev_0bk t(1)];
        t=textread('0bk_places.txt');
        ev_0bk=[ev_0bk t(1)];
        t=textread('2bk_body.txt');
        ev_2bk=[ev_2bk t(1)];
        t=textread('2bk_tools.txt');
        ev_2bk=[ev_2bk t(1)];
        t=textread('2bk_faces.txt');
        ev_2bk=[ev_2bk t(1)];
        t=textread('2bk_places.txt');
        ev_2bk=[ev_2bk t(1)];
        
        cd([basedir n1{pdx} '\EVs_RL'])
        
        t=textread('0bk_body.txt');
        ev_0bk=[ev_0bk t(1)];
        t=textread('0bk_tools.txt ');
        ev_0bk=[ev_0bk t(1)];
        t=textread('0bk_faces.txt');
        ev_0bk=[ev_0bk t(1)];
        t=textread('0bk_places.txt');
        ev_0bk=[ev_0bk t(1)];
        t=textread('2bk_body.txt');
        ev_2bk=[ev_2bk t(1)];
        t=textread('2bk_tools.txt');
        ev_2bk=[ev_2bk t(1)];
        t=textread('2bk_faces.txt');
        ev_2bk=[ev_2bk t(1)];
        t=textread('2bk_places.txt');
        ev_2bk=[ev_2bk t(1)];
        
        
        %CONVERT time to TR, TR=0.72 seconds. 
        ev_0bk=floor(ev_0bk/0.720);
        ev_2bk=floor(ev_2bk/0.720);
        ev_0bk(5:8) = ev_0bk(5:8)+405;
        ev_2bk(5:8) = ev_2bk(5:8)+405;
        ev_0bk=sort(ev_0bk);
        ev_2bk=sort(ev_2bk);
                
        onsets(1:8,1:2, pdx) = [ev_0bk' ev_2bk'];
    end
end



%%
for pdx = 1:length(n2)
    try
        cd(basedir)
        dat2(:,:,pdx)=csvread([n2{pdx} '_WM_merged_atlas_glasserROIs.csv']);
        
        % get even onset data
        cd([basedir n2{pdx} '\EVs_LR'])
        
        ev_0bk=[]; ev_2bk=[];
        t=textread('0bk_body.txt');
        ev_0bk=[ev_0bk t(1)];
        t=textread('0bk_tools.txt ');
        ev_0bk=[ev_0bk t(1)];
        t=textread('0bk_faces.txt');
        ev_0bk=[ev_0bk t(1)];
        t=textread('0bk_places.txt');
        ev_0bk=[ev_0bk t(1)];
        t=textread('2bk_body.txt');
        ev_2bk=[ev_2bk t(1)];
        t=textread('2bk_tools.txt');
        ev_2bk=[ev_2bk t(1)];
        t=textread('2bk_faces.txt');
        ev_2bk=[ev_2bk t(1)];
        t=textread('2bk_places.txt');
        ev_2bk=[ev_2bk t(1)];
        
        cd([basedir n2{pdx} '\EVs_RL'])
        
        t=textread('0bk_body.txt');
        ev_0bk=[ev_0bk t(1)];
        t=textread('0bk_tools.txt ');
        ev_0bk=[ev_0bk t(1)];
        t=textread('0bk_faces.txt');
        ev_0bk=[ev_0bk t(1)];
        t=textread('0bk_places.txt');
        ev_0bk=[ev_0bk t(1)];
        t=textread('2bk_body.txt');
        ev_2bk=[ev_2bk t(1)];
        t=textread('2bk_tools.txt');
        ev_2bk=[ev_2bk t(1)];
        t=textread('2bk_faces.txt');
        ev_2bk=[ev_2bk t(1)];
        t=textread('2bk_places.txt');
        ev_2bk=[ev_2bk t(1)];
                
        %CONVERT time to TR
        ev_0bk=floor(ev_0bk/0.720);
        ev_2bk=floor(ev_2bk/0.720);
        ev_0bk(5:8) = ev_0bk(5:8)+405;
        ev_2bk(5:8) = ev_2bk(5:8)+405;
        ev_0bk=sort(ev_0bk);
        ev_2bk=sort(ev_2bk);
        
        onsets2(1:8,1:2, pdx) = [ev_0bk' ev_2bk'];
    end
end

%%
% stack the events for all partiicpants in order to make averages. 
% stacking it is a bit hacky, but it'll do for a descriptive. Formal
% statistical analysis via FIR may be preferable, but is hard for me to
% impliment because of how the HCP GLms are set up and how the data is
% stored on our severers ( I can't write to the directories properly and
% copying all files is very space intensive)

% group 1, low activators
for rdx=1:360
    n=1;
   for pdx=1:100
        for idx=1:8
            try
                o=onsets(idx, 1,pdx);
                stacked_0bk(n, :,rdx)= squeeze(dat(rdx, o-2:o+dur,pdx ));
                
                o=onsets(idx, 2,pdx);
                stacked_2bk(n, :,rdx)= squeeze(dat(rdx, o-2:o+dur,pdx ));
                
                n=n+1;
            end
        end
    end
    n1_mean_2bk(:,rdx) = mean(stacked_2bk(:,:,rdx))';
    n1_mean_2bk(:,rdx) = n1_mean_2bk(:,rdx) - mean(n1_mean_2bk(1:5,rdx));
    
    n1_mean_0bk(:,rdx) = mean(stacked_0bk(:,:,rdx))';
    n1_mean_0bk(:,rdx) = n1_mean_0bk(:,rdx)  - mean(n1_mean_0bk(1:5,rdx));
end

% group 2
for rdx=1:360
    n=1;
   for pdx=1:100
        for idx=1:8
            try
                o=onsets2(idx, 1,pdx);
                stacked_0bk2(n, :,rdx)= squeeze(dat2(rdx, o-2:o+dur,pdx ));
                
                o=onsets2(idx, 2,pdx);
                stacked_2bk2(n, :,rdx)= squeeze(dat2(rdx, o-2:o+dur,pdx ));
                
                n=n+1;
            end
        end
    end
    n2_mean_2bk(:,rdx) = mean(stacked_2bk2(:,:,rdx))';
    n2_mean_2bk(:,rdx) = n2_mean_2bk(:,rdx) - mean(n2_mean_2bk(1:5,rdx));
    
    n2_mean_0bk(:,rdx) = mean(stacked_0bk2(:,:,rdx))';
    n2_mean_0bk(:,rdx) = n2_mean_0bk(:,rdx)  - mean(n2_mean_0bk(1:5,rdx));
 end




%%

rois=[263 83 63 111 181 209]

for r=rois
figure; plot([n1_mean_2bk(:,r) n1_mean_0bk(:,r)]); hold on
plot([n2_mean_2bk(:,r) n2_mean_0bk(:,r)], 'linewidth', 3)
end


%[n1_mean_2bk(:,r) n1_mean_0bk(:,r) n2_mean_2bk(:,r) n2_mean_0bk(:,r)]


