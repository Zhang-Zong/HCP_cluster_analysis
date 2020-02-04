function hcp_cognitive_anovas
% load clustering data (Wards) and cognitive data
% do one-way anovas for differences in cognition across clusters, for solutions
% ranging from 2 to 10
% 
% runs for cogntivie test scores, as well as task performance. All results FDR corrected.

load hcp_cog_dat
load hcp_hiearechical

% relationships between clusters (ward's, full sample) and cogntivie scores
for idx = 1:12
    for mdx = 1:6
        for cdx = 1:10
            C = cluster(Zall(:, :,mdx), 'maxclus', cdx);
            msize(mdx,cdx) = min(histc(C,1:cdx)); %number of people in the smallest cluster
            tinc(1:822) = 1;
            if msize(mdx,cdx) < 20  % exclude any clsuters with less than 20 members
                mn= msize(mdx,cdx);
                while mn < 20
                    ex = find(histc(C,1:max(C)) == mn);
                    tinc(squeeze(C == ex)) = 0;
                    C=C(tinc==1);
                    C(C> ex) = C(C>ex)-1;
                    mn =   min(histc(C,1:max(C)));
                end
            end
            [pval_wards(mdx,cdx,idx) atab(:,:,mdx,cdx,idx)] = anova1(zcog_822(tinc==1,idx),C, 'off');
        end
    end
end
% calc FDR and make anything which doesn't pass significance = 1
[p_fdr, p_masked] = fdr( nonzeros(pval_wards(:)), 0.05);
pval_wards(pval_wards > p_fdr) = 1; 


% relatipnships between clusters (ward's, full sample) and in scanner fmri 
% task performance metrics
for idx = 1:24
    for mdx = 1:6
        for cdx = 2:10
            C = cluster(Zall(:, :,mdx), 'maxclus', cdx);
            msize(mdx,cdx) = min(histc(C,1:cdx)); %number of people in the smallest cluster
            tinc(1:822) = 1;
            if msize(mdx,cdx) < 20  % exclude any clsuters with less than 20 members
               mn= msize(mdx,cdx);
               cc=C; 
                while mn < 20
                    ex = find(histc(C,1:max(C)) == mn);
                    tinc(squeeze(C == ex)) = 0;
                    cc(cc==ex) = 0;
                    mn =   min(nonzeros(histc(cc(cc>0),1:max(cc))));
                end
                C=C(tinc==1);
            end
            pval_tasks(mdx,cdx,idx) = anova1(ztask_822(tinc==1,idx),C, 'off');
        end
    end
end
% calc FDR and make anything which doesn't pass significance = 1
[p_fdr, p_masked] = fdr( nonzeros(pval_tasks(:)), 0.05);
pval_tasks(pval_tasks > p_fdr) = 1; 


save cognitive_anovas pval_wards pval_wards zcog_822 ztask_822


%% 
%by cogitive scores, but randomly permute the cluster membership, to
%convince me we don't find a bunch of false positives
% code fragment I kept for future use and/or interest
for idx = 1:12
    for mdx = 1:6
        for cdx = 2:10
            C = cluster(Zall(:, :,mdx), 'maxclus', cdx);
            C=C(randperm(length(C)));
            msize(mdx,cdx) = min(histc(C,1:cdx)); %number of people in the smallest cluster
            tinc(1:822) = 1;
            if msize(mdx,cdx) < 20  % exclude any clsuters with less than 20 members
                mn= msize(mdx,cdx);
                cc=C; 
                while mn < 20
                    ex = find(histc(C,1:max(C)) == mn);
                    tinc(squeeze(C == ex)) = 0;
                    cc(cc==ex) = 0;
                    mn =   min(nonzeros(histc(cc(cc>0),1:max(cc))));
                end
                C=C(tinc==1);
            end
            tpval(mdx,cdx,idx) = anova1(zcog_822(tinc==1,idx),C, 'off');
        end
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
t(:,1,:)= 0; %k=1, not possible, no result
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

