function hcp_group_clus(mod, name, C)
% hcp_group_clus(mod, name, C)
%
% Run a group analysis on HCP data for a given cluster solution, greating a
% group map for each cluster in that modality. 
%
% mod is a text of modality, name is a lsit of names, C is the cluster
% memberships


basedir  =  'F:\HCP900\data8\';
outdir = 'F:\HCP900\group8\all_wards\';

numclus=max(C);
clusname = ['wards_' mod '_clus' num2str(numclus)]
clus_var = C;

load group_t12.mat
con_names = {'c1'; };
confile = 'cope1.nii'

for clus = 1:max(clus_var)
    
    condir = [outdir clusname '_' num2str(clus) '/' ];
    mkdir(condir)
    job.dir = {[condir '/']};
    cd(condir);
    %clear scans in current batch
    job.des.t1.scans = {};
    n=1;
    for pdx = 1:length(clus_var)
        if clus_var(pdx) == clus
            job.des.t1.scans(n) = {[basedir name{pdx} '/' mod '/' confile  ',1']};
            n=n+1;
        end
    end
    
    try
        if ~isempty(ls('SPM.mat'))
            delete SPM.mat
        end;
    end
    
    job.masking.tm.tm_none=1;
    job.masking.im= 0;
    job.masking.em= {''};
    job.globalc.g_omit=1;
    job.globalm.glonorm=1;
    job.globalm.gmsca.gmsca_no=1;
    
    spm_run_factorial_design(job)
    
    load SPM
    SPM=spm_spm(SPM)
    save('SPM', 'SPM')
    
    cons = [1; -1];
    names = {'pos', 'neg'};
    curdir = pwd;
    analyze_spm_contrasts(curdir, cons, names);
    
end
