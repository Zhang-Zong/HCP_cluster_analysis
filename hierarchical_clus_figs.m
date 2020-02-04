function hierarchical_clus_figs

load hcp_data8.mat

cd F:\HCP900\figs\large_dendros

% make dendrogram. Not used in paper, not hyper informative, but nice to
% see. 
for mdx=1:6
    figure; set(gcf, 'Position', [-1097.4         93.8       3171.2          972]); 
    [h t pe] = dendrogram(Zall(:,:,mdx), 0, 'reorder', hord(:,mdx)); 
    set(gca,'XTick',[], 'YTick', []); 
    set( h, 'Color', 'k', 'Linewidth', 2 );
    saveas(gcf,[mods{mdx} '_HCwards_822.tiff'])
    
end

%%
% run group analysis code, used in Figure 1 and supplemental figs 1-8
% data is actually opend in wb_view, after converting the spmT files back
% into dscalar. 
for mdx=1:6
    for cdx=2:6
        C=cluster(Zall(:,:,mdx), 'Maxclus', cdx);
        hcp_group_clus(mods{6}, name_inc, C)
    end
end

%%
% FIGURE 2
% make color figs of PValues from this analysis, presenting (log-10)
% p-vales for significant relationships as colors (yellow-red)
%
% need to run hcp_cognitive_anovas.m first

clear out
t=log10(pval_wards)*-1;
t(:,1,:)= 0;
for cdx = 2:10
    n=1
    for idx = 1:12
        for  mdx=1:6
            %out(idx,cdx,mdx) = pval_wards(mdx,cdx,idx);
            out(idx,cdx,mdx) = t(mdx,cdx,idx)
            n=n+1;
        end
    end
end

%actual images used in Figure 2
cd F:\HCP900\figs\cog
%make non significant a negative number to increase visibility on figure,
%make sure it is clear what is and is not significant
out(out==0) = -5;
% hot color map excluding black for prettyness
h=hot; 
hmap = (h(end:-1:20,:));
for mdx=1:6
    
    figure;  imagesc(out(:,2:10,mdx), [-2 10]);
    colormap(hmap); 
    set(gcf, 'Position', [768.2 449 265.6 420])
    set(gca,'XTick',0.5:8.5, 'YTick', 0.5:1:11.5)
    grid on
    set(gca,'XTicklabel',[], 'YTicklabel', []);
    saveas(gcf,[mods{mdx} '_zcog_logp_822.tiff'])
end

        
%% 
% agrement matricies, hierarchical
% FIGURE 3 data
load wards_agree.mat  

% In order to make sense of the clsuter probability matiecise (here called
% 'agree', ebcaue i origonally calle dthem agrement matricies), we need to
% run optimalleaforder. This sorts the matrix in a way that makes visual
% sense and is infmative. 

for mdx=1:6 % modality/scan
    for cdx = 2:10
        [mdx cdx]
        a=agree(:,:,cdx,mdx);
        ad=pdist(a);
        Z = linkage(a, 'average');
        word(1:822,cdx,mdx) = optimalleaforder(Z,ad,'criteria', 'group');
    end
end

% plot and save cluste prob matricies, data used in FIGURE 3
cd F:\HCP900\figs\agree
for mdx=1:6
    for cdx = 2:10  
        a=agree(:,:,cdx,mdx);
        figure; imagesc(a(word(:,cdx,mdx),word(:,cdx,mdx)), [0 1]);
        colormap(hot)
        
        set(gca,'XTicklabel',[], 'YTicklabel', []);
        
        saveas(gcf,['zzHC_' mods{mdx} '_' num2str(cdx) '.tif'])
    end
end

%%
%example agree matrix components, used in Figure 4

mdx=6
cdx=4
a=agree(:,:,cdx,mdx);
figure; imagesc(a(ord(:,cdx,mdx),ord(:,cdx,mdx)), [0 1]);
colormap(hot)
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] =  pca(agree(:,:,cdx,mdx)); %pca(sdis(:,:,mdx));
figure; imagesc(SCORE(ord(:,cdx,mdx),1));colormap(hot)
figure; imagesc(SCORE(ord(:,cdx,mdx),2));colormap(hot)
figure; imagesc(SCORE(ord(:,cdx,mdx),3));colormap(hot)

% component variance explained in figure 4 are from exp below, I made the
% plot in excel, super obvously. 
for mdx=1:6
    for cdx=2:10
        [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] =  pca(agree(:,:,cdx,mdx));
        agree_comps(:,1:3,mdx,cdx)=SCORE(:,1:3) ;
        exp(:, mdx,cdx)=EXPLAINED;
    end
end

%%
% FIGURE 5, this is where shit gets super cool. 


load wards_agree.mat 

% first lets plot the component spaces, the 'snake' and 'tortilla' plot. My
% favoites! This is for the 4 clsuter solution presented in Figure 5. 

% clusters order as in group fig; lowest/negative to highest/pos. The order
% of clsuters in hierarcical is arbitrary so we need to reorganize now for
% consistency. 
% order in group figure
cluso = [
    4 3 1 2
    2 1 4 3
    4 2 3 1
    3 4 2 1
    1 2 3 4
    2 1 4 3];
     
clr =  'kbgr'
for mdx=1:6
    cdx=4
    gh = figure; hold on; grid on; 
        
    [COEFF, SCORE, latent, tsq, explain] =  pca(agree(:,:,cdx,mdx)); 
    ss=SCORE(:,1:3) ;
    c = cluster(Zall(:,:,mdx),'MaxClust',cdx);
    for idx = 1:max(c);
         scatter3(ss(c==cluso(mdx,idx),1), ss(c==cluso(mdx,idx),2), ss(c==cluso(mdx,idx),3), 'o', clr(idx))
    end
    set(gca,'XTickLabel',[], 'YTickLabel', [], 'ZTickLabel', []);
    ax=gca; ax.GridAlpha = 1;
    rotate3d ON
end
% I did not auto save, I manually rotated each figure to a position I liked
% and saved it that way. 

% now we need the regression maps, Figure 5, central panel. For this we
% need to run an SPM group analysis, and then convert the outputs back to
% dscalars, and visualized using wb_view. 

% this one script generates all the maps we need. 
hcp_group_clus_regress


%% 

% FIGURE 6 radar plots by network

% Unfortunatly I didn't retain this peice of code. I extracted the numbers
% using Erin Dickie's ciftify toolbox (https://github.com/edickie/ciftify)
% but using the yeo 7 network map as an atlas, and looping through t-maps
% with ciftify-meants. 
%
% graphs were made in excel. 


%%

% Correlations with cognitive scores and cluster bootstrap components. 
% FIGURE 7! 
% agree_comps calcualted above, need to run this. 

% corr between components and cog scores, using k=4. 
for mdx=1:6
    cdx=4
    [rr pp] = corr(agree_comps(:,1:3,mdx,cdx), cog_822, 'rows', 'pairwise', 'type', 'Spearman')
    ss_cor(:,:,mdx) = rr; 
    cc_corp(:,:,mdx) = pp; 
end

% fdr correction, change non sig correlations to 0. 
p_fdr= fdr(cc_corp(:), 0.05)
ss_cor(cc_corp > p_fdr) = 0

% since the components have no direction per se, we will express all
% correlations as absolute values
ss_cor=abs(ss_cor);

%plot correlations as heat maps
cd F:\HCP900\figs\cog
h=hot; 
% reverse the ehat map to make this scaled as I wish. 
hmap = (h(end:-1:20,:));

for mdx=1:6
    figure;  imagesc(ss_cor(:,:,mdx)', [0 0.25]);
    colormap(hmap)
    set(gcf, 'Position', [768.2 449 265.6 420])
    set(gca,'XTick',0.5:2.5, 'YTick', 0.5:1:11.5)
    grid on
    set(gca,'XTicklabel',[], 'YTicklabel', []);
    saveas(gcf,['agreeSS_' mods{mdx} '_cogcorrs.tiff'])
end

%% 
% run psuedo-simulations as in Figure 8. 

for cdx=4 % all sims on k=4 data
    for mdx=1:6
        
        demean_dat = detrend(data(:,:,mdx)','constant')';
        
        % simulation of demeaned data
        [ARI, demean_agree(:,:,cdx,mdx), CRI, Cout] = cluster_bootstrap(demean_dat, cdx);
        
        %simulate shuffled data, not demeaned
        shuff_dat = [];
        for pdx=1:822
            shuff=nonzeros(reshape(sb(randperm(590),:)', 1, 59000)); % random shuffle in blocks of 100 verticies
            shuff_dat(pdx,:) = data(pdx, shuff,mdx);
        end
        [ARI, shuff_agree(:,:,cdx,mdx), CRI, Cout] = cluster_bootstrap(shuff_dat, cdx);
        % demeaned AND shuffled!
        demshuff_dat = detrend(shuff_dat','constant')';
        [ARI, demshuff_agree(:,:,cdx,mdx), CRI, Cout] = cluster_bootstrap(demshuff_dat, cdx)  ;
    end
end

% this is the last one, randomized with simulated networks added. 
cdx=4
for mdx=1:5
    
        shuff_dat = [];
        for pdx=1:822
            shuff_dat(pdx,:) = data(pdx, randperm(size(data,2)),mdx);
        end
        demshuff_dat = detrend(shuff_dat','constant')';
        % the simulated networks, n1m n2, and n3
        n1 = 5000:16000; n2 = 25000:31000; n3 = 40000:46000; 
        
        rat = (1/822:1/822:1)+1; %amout of network activity, from lower to higher, to simulate modulating networks
        %random order, to apply to network modualtions, so each of the
        %three is independantly modualted
        ro = [randperm(822)' randperm(822)' randperm(822)']; 
        
        % add the sim networks to the demenaned and shuffled data
        for pdx=1:822 
            demshuff_dat(pdx, n1) = abs(demshuff_dat(pdx, n1)) * rat(ro(pdx,1));
            demshuff_dat(pdx, n2) = abs(demshuff_dat(pdx, n2)) * rat(ro(pdx,2))*-1;
            demshuff_dat(pdx, n3) = abs(demshuff_dat(pdx, n3)) * (rat(ro(pdx,3))-1.5)*2;
        end
       
        [ARI, sim_mean_n3_agree(:,:,cdx,mdx), CRI, Cout] = cluster_bootstrap(demshuff_dat, cdx, 500, .75)  ;
     
end

% MAKE PLOTS FOR FIGURE 8!!
cdx=4;
% top row, demeaned
for mdx=1:6
    aa=demean_agree(:,:,cdx,mdx);
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] =  pca(aa);
    figure; scatter3(SCORE(:,1), SCORE(:,2),SCORE(:,3), 10, 'k'); 
    set(gca,'XTicklabel',[], 'YTicklabel', [], 'ZTicklabel', []);
end
% second row, shuffled, not demeaned
for mdx=1:6
    aa=shuff_agree(:,:,cdx,mdx);
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] =  pca(aa);
    figure; scatter3(SCORE(:,1), SCORE(:,2),SCORE(:,3), 10, 'k'); 
    set(gca,'XTicklabel',[], 'YTicklabel', [], 'ZTicklabel', []);
end
% third row, random noise, shuffled and demeaned
for mdx=1:6
    aa=demshuff_agree(:,:,cdx,mdx);
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] =  pca(aa);
    figure; scatter3(SCORE(:,1), SCORE(:,2),SCORE(:,3), 10, 'k'); 
    set(gca,'XTicklabel',[], 'YTicklabel', [], 'ZTicklabel', []);
end
% simulated networks, fourth rown, figure 8
for mdx=1:6
    aa=sim_mean_n3_agree(:,:,cdx,mdx);
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] =  pca(aa);
    figure; scatter3(SCORE(:,1), SCORE(:,2),SCORE(:,3), 10, 'k'); 
    set(gca,'XTicklabel',[], 'YTicklabel', [], 'ZTicklabel', []);
end

%%
 % Data from figure 9 is made with a separate function, 
WM_SS1_hrf_EXTRACTIONS

% the code is a mess, but it worked. NO APPLOGIES! 


%% SUPPLEMENTAL FIGURE ON PMAT BOXPLOTS
% SUP FIG 11

cluso = [ % order of clusters, neg to pos
    4 3 1 2
    2 1 4 3
    4 2 3 1
    3 4 2 1
    1 2 3 4
    2 1 4 3];

for mdx=1:6
    cdx=4
    clear CC
    c = cluster(Zall(:,:,mdx),'MaxClust',cdx)*10;
    for idx=1:4
        CC(c==idx*10) = cluso(mdx,idx);
    end
    anova1(zcog_822(:,5),CC')
    title(mods{mdx})
    set(gcf, 'Position', [144.33 332.33 330.67 399.33])
    set(findobj(gca,'type','line'),'linew',2)
end

%%
% SUP FIG 12
clr='rgcbmkrgcbmk'
cd F:\HCP900\figs/agree_comps
for mdx=1:6
    cdx=2:10
    figure; hold on; grid on;
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] =  pca(agree(:,:,cdx,mdx));
    ss=SCORE(:,1:3) ;
    c = cluster(Zall(:,:,mdx),'MaxClust',cdx);
    for idx = 1:max(c);
        scatter3(ss(c==idx,1), ss(c==idx,2), ss(c==idx,3), 'o', clr(idx))
    end
    view(50,20)
    %saveas(gcf,['agreeSS_' mods{mdx} '_c' num2str(cdx) '_3Dcomps.tiff'])
end


%%
% SUP FIG 13
%% plot colors by mean activity for tortilla/snakes
clr=parula(822); % color map by subs
for mdx=1:6
    cdx=4;
    
    mm = mean(data(:,:,mdx),2); 
    sd = std(data(:,:,mdx)')';
    
    CV = sd ./ mm  ;
    
    ss=sortrows([mm (1:822)']);
    ss=ss(:,2); % data ordered according to mean tstat across the brain 
    
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] =  pca(agree(:,:,cdx,mdx));
    figure; scatter3(SCORE(ss,1), SCORE(ss,2),SCORE(ss,3), 10, clr)
    
end

%%
% SUP FIG 14
% multidimentional scaling

for mdx=1:6
    mdx
    cdx=4
    D=1-agree(:,:,cdx,mdx);
    MD(1:822, 1:10, mdx) = mdscale(D,10); 
    figure; scatter3(MD(:, 1, mdx), MD(:, 2, mdx),MD(:, 3, mdx), 10, 'k');
end

%%
% SUP FIG 15
% HISTORGRAMS, supplemental figure 15
for mdx=1:6
    figure; hist(agree_comps(:,1,mdx,4), 40)
end

%%
% SUP FIG 16
% Corrleate component scores with fMRI task scores

%%%%%%%%%%%%%%%%%%%%%%
% fmri Task data
for mdx=1:6
    cdx=4
    [rr pp] = corr(agree_comps(:,1:3,mdx,cdx), task_822, 'rows', 'pairwise', 'type', 'Spearman');
    ss_task(:,:,mdx) = rr; 
    ts_corp(:,:,mdx) = pp; 
end

% fdr correction, change non sig correlations to 0. 
p_fdr= fdr(ts_corp(:), 0.05)
ss_task(ts_corp > p_fdr) = 0

% since the components have no direction per se, we will express all
% correlations as absolute values
ss_task=abs(ss_task);

%plot correlations as heat maps
cd F:\HCP900\figs\cog
h=hot; 
hmap = (h(end:-1:20,:));
for mdx=1:6
    
    figure;  imagesc(ss_task(:,:,mdx)', [0 0.25]);
    colormap(hmap)
    set(gcf, 'Position', [768.2 449 265.6 420])
    set(gca,'XTick',0.5:2.5, 'YTick', 0.5:1:23.5)
    grid on
    set(gca,'XTicklabel',[], 'YTicklabel', []);
   saveas(gcf,['agreeSS_' mods{mdx} '_taskcorrs.tiff'])
end

































%%




%%


% clusters order as in group fig; lowest/negative to highest/pos

%order in group figure
cluso = [
    4 3 1 2
    2 1 4 3
    4 2 3 1
    3 4 2 1
    1 2 3 4
    2 1 4 3];
     
clr =  'kbgr'
for mdx=1:6
    cdx=4
    gh = figure; hold on; grid on; 
        
    [COEFF, SCORE, latent, tsq, explain] =  pca(agree(:,:,cdx,mdx)); 
    ss=SCORE(:,1:3) ;
    c = cluster(Zall(:,:,mdx),'MaxClust',cdx);
    for idx = 1:max(c);
         scatter3(ss(c==cluso(mdx,idx),1), ss(c==cluso(mdx,idx),2), ss(c==cluso(mdx,idx),3), 'o', clr(idx))
    end
    set(gca,'XTickLabel',[], 'YTickLabel', [], 'ZTickLabel', []);
    ax=gca; ax.GridAlpha = 1;
    rotate3d ON
end

%% plot colrs for cognitive ability

ct = cog_822(:,5)
ct(isnan(ct)==1) = floor(mean(ct,'omitnan'));
cmap= parula(max(ct)-min(ct)+1);

for mdx=1:6
    cdx=4
    gh = figure; hold on; grid on; 
        
    [COEFF, SCORE, latent, tsq, explain] =  pca(agree(:,:,cdx,mdx)); 
    ss=SCORE(:,1:3) ;
    c = cluster(Zall(:,:,mdx),'MaxClust',cdx);
    for idx = 1:822
        tcc(idx, 1:3) = cmap(ct(idx) - min(ct)+1,:);
         %scatter3(ss(idx,1), ss(idx,2), ss(idx,3), cmap(ct(idx) - min(ct)+1,:))
    end
    scatter3(ss(:,1), ss(:,2), ss(:,3), ones(822,1)*10, tcc)
    set(gca,'XTickLabel',[], 'YTickLabel', [], 'ZTickLabel', []);
    ax=gca; ax.GridAlpha = 1;
    rotate3d ON
end


scatter3(ss(:,1), ss(:,2), ss(:,3),  (ct-min(ct)), 'k')



%% plot colors by mean activity for tortilla/snakes


clr=parula(822); % color map by subs
for mdx=1:6
    cdx=4;
    
    mm = mean(data(:,:,mdx),2); 
    sd = std(data(:,:,mdx)')';
    
    CV = sd ./ mm  ;
    
    ss=sortrows([mm (1:822)']);
    ss=ss(:,2); % data ordered according to mean tstat across the brain 
    
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] =  pca(agree(:,:,cdx,mdx));
    figure; scatter3(SCORE(ss,1), SCORE(ss,2),SCORE(ss,3), 10, clr)
    
end





%%
% clusters plotted along first 2 components of PCA of raw data matrix
% not included in paper, makes data clouds. Interesting but in the end
% rejected as not a terribly useful way to visualize the data. 
cd F:\HCP900\figs\pca_plots
for mdx=1:6
    clr='bgkrm' % color order
    shp='o*svhx+' %shape order
    nclus=4; %number of clsuters to plot
   [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] =  pca(data(:,:,mdx)); 
   c  = cluster(Zall(:,:,mdx), 'MaxClust',nclus);
   figure; hold on
   for cdx=1:nclus
       plot(SCORE(c==cdx, 1), SCORE(c==cdx, 2), [clr(cdx) shp(cdx)]);
   end
   set(gca,'XTick',[], 'YTick', []);
   saveas(gcf,[mods{mdx} '_wards_pcaplots_822.tiff'])
end
       




