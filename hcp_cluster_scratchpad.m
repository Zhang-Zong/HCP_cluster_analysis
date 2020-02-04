% this is a scratchpad script which contains many of the calls and commands
% used to clsuter the HCP data


%%
% Basic clustering variables and parameters
Ze= linkage(data(:,:,1), 'ward');
Zg= linkage(data(:,:,2), 'ward');
Zl= linkage(data(:,:,3), 'ward');
Zr= linkage(data(:,:,4), 'ward');
Zs= linkage(data(:,:,3), 'ward');
Zw= linkage(data(:,:,6), 'ward');

figure; [h t pe] = dendrogram(Ze, 0,'ColorThreshold', 2000);
figure; [h t pg] = dendrogram(Zg, 0,'ColorThreshold', 2000);
figure; [h t pl] = dendrogram(Zl, 0,'ColorThreshold', 2000);
figure; [h t pr] = dendrogram(Zr, 0,'ColorThreshold', 2000);
figure; [h t ps] = dendrogram(Zs, 0,'ColorThreshold', 2000);
figure; [h t pw] = dendrogram(Zw, 0,'ColorThreshold', 2000);

%very provisonal number of clsuters based on visiual examination of the data
Ce = cluster(Ze,'MaxClust',6);
Cg = cluster(Zg,'MaxClust',5);
Cl = cluster(Zl,'MaxClust',2);
Cr = cluster(Zr,'MaxClust',3);
Cs = cluster(Zs,'MaxClust',3);
Cw = cluster(Zw,'MaxClust',5);

Call = [Ce Cg Cl Cr Cs Cw];

Ce = cluster(Ze,'MaxClust',3);
Cg = cluster(Zg,'MaxClust',3);
Cl = cluster(Zl,'MaxClust',3);
Cr = cluster(Zr,'MaxClust',3);
Cs = cluster(Zs,'MaxClust',3);
Cw = cluster(Zw,'MaxClust',3);

Call3 = [Ce Cg Cl Cr Cs Cw];

Ce = cluster(Ze,'MaxClust',2);
Cg = cluster(Zg,'MaxClust',2);
Cl = cluster(Zl,'MaxClust',2);
Cr = cluster(Zr,'MaxClust',2);
Cs = cluster(Zs,'MaxClust',2);
Cw = cluster(Zw,'MaxClust',2);

Call2 = [Ce Cg Cl Cr Cs Cw];

Zall(:,:,1)=Ze; 
Zall(:,:,2)=Zg;
Zall(:,:,3)=Zl;
Zall(:,:,4)=Zr;
Zall(:,:,5)=Zs;
Zall(:,:,6)=Zw;

pall = [ pe;   pg ;  pl;   pr;   ps;   pw ];

for mdx = 1:6
    mdx
    pdis = pdist(data(:,:,mdx), 'euclidean');
    sdis(:,:,mdx) = squareform(pdis); 
end


for mdx = 1:6
    m1 = min(nonzeros(sdis(pall(mdx,:), pall(mdx,:),mdx))); 
    m2 = max(nonzeros(sdis(pall(mdx,:), pall(mdx,:),mdx))); 
    figure; imagesc(sdis(pall(mdx,:), pall(mdx,:),mdx), [m1*1.2 m2*.7])
end


%%
% ARI for main cluster solutions
for cdx = 2:10
    for mdx = 1:6
        C = cluster(Zall(:,:,mdx),'MaxClust',cdx);
        for rdx = 1:6
            C2 = cluster(Zall(:,:,rdx),'MaxClust',cdx);
            [AR(mdx,rdx),RI,MI,HI]=RandIndex(C, C2)
        end
    end
end

%%
% Enter-fix behav data
for pdx=1:length(name_inc)
    id = str2num(name_inc{pdx});
%     em_rt(pdx) = emotion_rt(find(cname == id));
%     er40_rt(pdx) = t(find(cname == id));
    er40_cr(pdx) = tt(find(cname == id));
end

%%
% ANOVA of differences in cog between clsuters
for mdx = 1:6
    for cdx = 2:10
        C = cluster(Zall(:,:,mdx),'MaxClust',cdx);
        tinc(1:822) = 1;
        msize(mdx,cdx) = min(histc(C,1:cdx)); %number of people in the smallest cluster
        if msize(mdx,cdx) < 20  % exclude any clsuters with less than 20 members
            mn= msize(mdx,cdx);
            while mn < 20
                ex = find(histc(C,1:max(C)) == mn);
                tinc(squeeze(C == ex)) = 0;
                C=C(tinc==1);
                C(C> ex) = C(C>ex)-1;
                mn = min(histc(C,1:max(C)));
            end
        end
        pval(mdx,cdx,1) = anova1(cog_pca2(tinc==1,1), C, 'off');
        pval(mdx,cdx,2) = anova1(cog_pca2(tinc==1,2), C, 'off');
        pv_attn(mdx,cdx) = anova1(cog_all(tinc==1,9), C, 'off');
    end
end
msize
pval


%%

% group maps for all clsuter solutions

for mdx=1:6
    for cdx=2:10
        C=cluster(Zall(:,:,mdx), 'maxClus', cdx)
         hcp_group_clus(mods{mdx}, name_inc, C)
    end
end

%%
% bootstrap clsuters

for mdx = 6
    for cdx = 5
        [ARI, perct_agree, CRI, Co] = cluster_bootstrap(data(:,:,6), cdx , 1000)
    end
end


%%
% sort agreement matricies
for mdx=1:6 % modality/scan
    for cdx = 2:10
        [mdx cdx]
        a=agree(:,:,cdx,mdx);
        ad=pdist(a);
        Z = linkage(a, 'average');
        ord(1:822,cdx,mdx) = optimalleaforder(Z,ad,'criteria', 'group');
    end
end

for cdx = 2:5
    a=agree(:,:,cdx,mdx);
    figure; imagesc(a(ord(:,cdx,mdx),ord(:,cdx,mdx)))
end


for mdx = 6
    figure; 
    for cdx = 2:10
        subplot(2,5,cdx)
          a=agree(:,:,cdx,mdx);
          imagesc(a(ord(:,cdx,mdx),ord(:,cdx,mdx)))
    end
end

for mdx=1:6
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(sdis(:,:,mdx));
ss=SCORE(:,1:3);
figure; scatter3(ss(:,1), ss(:,2), ss(:,3)); grid on
end

corr(cog_pca2, ss, 'rows', 'pairwise')

Y = quantile(cog(:,1),3)
cogq=ones(1,822)*4
for idx=3:-1:1
    cogq(cog(:,1) < Y(idx)) = idx;
end

clr='brgcmky'
figure; hold on;
for idx=1:4
     scatter3(ss(cogq==idx,1), ss(cogq==idx,2), ss(cogq==idx,3), 'o', clr(idx))
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 3D plots of agreement matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clr =  'rgbkcmyrgbk'
for mdx=1:6
    cdx=4
    
    gh = figure; hold on; grid on;
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] =  pca(agree(:,:,cdx,mdx)); 
    ss=SCORE(:,1:3) ;
    c = cluster(Zall(:,:,mdx),'MaxClust',cdx);
    for idx = 1:max(c);
         scatter3(ss(c==idx,1), ss(c==idx,2), ss(c==idx,3), 'o', clr(idx))
    end
    rotate3d ON
    
end


figure; hold
for idx = 1:max(c);
    plot(ss(c==idx,1), ss(c==idx,2), ['o' clr(idx)])
end

%%



for mdx = 1:6
    figure; plot(squeeze(mean(mean(ari2(:,:,2:10,mdx)))))
end


%% 
% cognitive data
n=1; 
for pdx = 1:899
    if inc(pdx) ==1
        cog_822(n,:) = cog_all(cogid == str2num(name{pdx}), :);
        n=n+1;
    end
end


%%

% multidimentional scaling


for mdx=1:6
    mdx
    cdx=4
    D=1-agree(:,:,cdx,mdx);
    MD(1:822, 1:3, mdx) = mdscale(D,3);    
end

mdx=2
Y=MD(:,1:3,mdx)
c=cluster(Zall(:,:,mdx), 'maxClus', cdx);
figure; hold
for idx = 1:4
    scatter3(Y(c==idx,1), Y(c==idx,2), Y(c==idx,3), 'o', clr(idx))
end



%%
cordis=pdist(data(:,:,6), 'euclidean');
scordis = squareform(cordis);

for idx = 1:822
    mcdis(idx) = mean(nonzeros(scordis(idx,:)));
end


s1=sortrows([ mcdis', (1:822)'], 'descend');

ord=[s1(:,2)]
figure; imagesc(scordis(ord,ord), [0.4 1]); colormap(winter)

figure; imagesc(scordis(ord,ord), [600 1400]); colormap(winter)


%%
% networks activity in SS regressions
cd F:\HCP900\group8\SS_yeo7_CSVs
for mdx=1:6
    for idx=1:4
        
        fname = ['HC_regress_SS_' mods{mdx} '_k4_c' num2str(idx) '_fsaverage.Yeo2011_7Networks_N1000.32k_fs_LR_meants.csv'];
        
        SS_nets(mdx, idx, :) = textread(fname)
        
    end
end
     

%%
% tacoplots colored by compoment 1
clr=hot(822)
for mdx=6
    cdx=4
    figure; hold on; grid on;
    ss=agree_comps(:,1:3, mdx, cdx);
    % reorder colormap by st component
    oo=sortrows([ss(:,1) (1:822)']); 
    oo=sortrows([oo(:,2), (1:822)'])
    scatter3(ss(:,1), ss(:,2), ss(:,3), ones(822,1)*10, clr(oo(:,2),:))
end

        