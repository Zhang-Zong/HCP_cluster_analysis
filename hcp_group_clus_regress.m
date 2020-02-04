% runs regression analyses of clsuter probability components, as in Figure 5, centre panel 
% we plot for all values of k, and I looked at them, but we only used k=4
% in the paper. 
%
% similar across other k's, but it gets messy complicated with too much
% information to handle. 

load wards_agree.mat 

for mdx=1:6
    
    for cdx=2:10
        [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] =  pca(agree(:,:,cdx,mdx)); %pca(sdis(:,:,mdx));
        ss=SCORE(:,1:3) ;
        
        basedir  =  'F:\HCP900\data8\';
        outdir = 'F:\HCP900\group8\';
        
        clusname = ['HC_regress_SS_' mods{mdx} '_k' num2str(cdx)]
        
        load group_t12.mat
        
        % set up regressors
        job.des = rmfield(job.des, 't1');
        job.des.mreg.incint = 1;
        
        job.des.mreg.mcov(1).cname = 'regressor1';
        job.des.mreg.mcov(1).iCC =1;
        
        job.des.mreg.mcov(2).cname = 'regressor2';
        job.des.mreg.mcov(2).iCC =1;
        
        job.des.mreg.mcov(3).cname = 'regressor3';
        job.des.mreg.mcov(3).iCC =1;
        
        job.des.mreg.mcov(1).c = ss(:,1);
        job.des.mreg.mcov(2).c = ss(:,2);
        job.des.mreg.mcov(3).c = ss(:,3);
        %regressors set
        
        con_names = {'c1'; };
        confile = 'cope1.nii'
        
        for clus = 1% jsut a residual from old code Iw as to lazy to remove. 
            
            condir = [outdir clusname]; % '_' num2str(clus) '/' ];
            mkdir(condir)
            job.dir = {[condir '/']};
            cd(condir);
            %clear scans in current batch
            job.des.mreg.scans = {};
            for pdx = 1:822
                job.des.mreg.scans(pdx) = {[basedir name_inc{pdx} '/' mods{mdx} '/' confile  ',1']};
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
            
            cons = [
                0 1 0 0;
                0 0 1  0;
                0 0 0 1;
                0 -1 0 0;
                0 0 -1 0;
                0 0 0 -1;];
            
            names = {'pos1', 'pos2';
                'pos3', 'neg1';
                'neg2', 'neg3';};
            
            curdir = pwd;
            analyze_spm_contrasts(curdir, cons, names);
            
        end
    end
end
