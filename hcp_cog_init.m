function hcp_cog_init
% a script for loading the cognitive data from the HCP study. Looks in the
% csv file HCP_subject_data.csv which includes relevant cognitive scores,
%  and task fMRI scores. Because of my understanding of the HCP data
%  sharing agreement, I have no shared this file :(

cd C:\Users\colin_hawco\Documents\GitHub\HCP_cluster_analysis
csvfname ='HCP_subject_data.csv'

load hcp_data8

cog_name = csvread(csvfname,1,0, [1  0 970 0]);
cog_dat = csvread(csvfname,1,2, [1  2 970 13]);
task_dat = csvread(csvfname,1,48, [1  48 970 71]);


% the cog csv doesn't line up with the imaging data, so we need to loop
% through names and find the matching index in the cog data
n=1; 
for pdx=1:length(name_inc)
    nm=str2num(name_inc{pdx}); 
    cog_822(n,1:12) = cog_dat(cog_name == nm, :);
    task_822(n,1:24) = task_dat(cog_name == nm, :);
    n=n+1; 
end

% z transform cog scores
for idx=1:12
    mn=mean(cog_822(:,idx), 'omitnan'); 
    sd=std(cog_822(:,idx), 'omitnan'); 
    zcog_822(:,idx) = (cog_822(:,idx) - mn) / sd; 
end
% z transform task scores
for idx=1:24
  mn=mean(task_822(:,idx), 'omitnan'); 
    sd=std(task_822(:,idx), 'omitnan'); 
    ztask_822(:,idx) = (task_822(:,idx) - mn) / sd; 
end  

% save variable names for later refrence
cog_var_names = {
    'PicSeq_AgeAdj'
    'CardSort_AgeAdj'
    'Flanker_AgeAdj'
    'ListSort_AgeAdj'
    'PMAT24_A_CR'
    'ReadEng_AgeAdj'
    'PicVocab_AgeAdj'
    'ProcSpeed_AgeAdj'
    'VSPLOT_CRTE'
    'SCPT_TP'
    'IWRD_TOT'
    'ER40_CRT'};

task_var_names = {
'Emotion_Task_Face_Acc'
'Emotion_Task_Shape_Acc'
'Emotion_Task_Face_Median_RT'
'Emotion_Task_Shape_Median_RT'
'Gambling_Task_Median_RT_Larger'
'Gambling_Task_Median_RT_Smaller'
'Language_Task_Story_Acc'
'Language_Task_Math_Acc'
'Language_Task_Story_Median_RT'
'Language_Task_Math_Median_RT'
'Language_Task_Story_Avg_Difficulty_Level'
'Language_Task_Math_Avg_Difficulty_Level'
'Relational_Task_Match_Acc'
'Relational_Task_Rel_Acc'
'Relational_Task_Match_Median_RT'
'Relational_Task_Rel_Median_RT'
'Social_Task_Random_Perc_Random'
'Social_Task_TOM_Perc_TOM'
'Social_Task_Random_Median_RT_Random'
'Social_Task_TOM_Median_RT_TOM'
'WM_Task_2bk_Acc'
'WM_Task_0bk_Acc'
'WM_Task_2bk_Median_RT'
'WM_Task_0bk_Median_RT' };

clear inc scvfname data idx inc mn sd mods n name name_inc nm pdx 
save hcp_cog_dat







