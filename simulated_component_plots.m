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


% MAKE PLTO FOR FIG 8
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







cdx=4
for mdx=1:5
    
        shuff_dat = [];
        for pdx=1:822
            shuff_dat(pdx,:) = data(pdx, randperm(size(data,2)),mdx);
        end
        
        demshuff_dat = detrend(shuff_dat','constant')';
        
        n1 = 5000:16000; n2 = 25000:31000; n3 = 40000:46000; 
        
        rat = (1/822:1/822:1)+1; %amout of network activity, from lower to higher
        ro = [randperm(822)' randperm(822)' randperm(822)'];
        
        for pdx=1:822
            demshuff_dat(pdx, n1) = abs(demshuff_dat(pdx, n1)) * rat(ro(pdx,1));
            demshuff_dat(pdx, n2) = abs(demshuff_dat(pdx, n2)) * rat(ro(pdx,2))*-1;
            demshuff_dat(pdx, n3) = abs(demshuff_dat(pdx, n3)) * (rat(ro(pdx,3))-1.5)*2;
        end
       
        [ARI, sim_mean_n3_agree(:,:,cdx,mdx), CRI, Cout] = cluster_bootstrap(demshuff_dat, cdx, 500, .75)  ;
     
end


for mdx=1:6
[COEFF, SCORE] =  pca(sim_mean_n3_agree(:,:,cdx,mdx));
figure; scatter3(SCORE(:,1), SCORE(:,2),SCORE(:,3), 10, 'k');  rotate3d on
set(gca,'XTicklabel',[], 'YTicklabel', [], 'ZTicklabel', []);
end

%% 


cdx=4
for mdx=6
    
        shuff_dat = [];
        for pdx=1:822
            shuff_dat(pdx,:) = data(pdx, randperm(size(data,2)),mdx);
        end
        
        demshuff_dat = detrend(shuff_dat','constant')';
        
        n1 = 5000:11000; n2 = 25000:31000; n3 = 40000:46000; 
        
        rat = (1/822:1/822:1)+1; %amout of network activity, from lower to higher
        ro = [randperm(822)' randperm(822)' randperm(822)'];
        
        for pdx=1:822
            demshuff_dat(pdx, n1) = abs(demshuff_dat(pdx, n1)) * rat(ro(pdx,1));
            demshuff_dat(pdx, n2) = abs(demshuff_dat(pdx, n2)) * rat(ro(pdx,2));
            demshuff_dat(pdx, n3) = abs(demshuff_dat(pdx, n3)) * rat(ro(pdx,3));
        end
        
        [ARI, sim_mean_n3pos_agree(:,:,cdx,mdx), CRI, Cout] = cluster_bootstrap(demshuff_dat, cdx, 500, .75)  ;
     
        
end


for mdx=1:6
[COEFF, SCORE] =  pca(sim_mean_n3_agree(:,:,cdx,mdx));
figure; scatter3(SCORE(:,1), SCORE(:,2),SCORE(:,3), 10, 'k');  rotate3d on
set(gca,'XTicklabel',[], 'YTicklabel', [], 'ZTicklabel', []);
end


% pure simulation using mu and std from WM
for pdx=1:822
    rd(pdx,:) = normrnd( 0.4721,0.4721,[1 size(data,2)]);
end

n1 = 5000:11000; n2 = 25000:31000; n3 = 40000:46000;

rat = (1/822:1/822:1)+1; %amout of network activity, from lower to higher
ro = [randperm(822)' randperm(822)' randperm(822)'];

for pdx=1:822
    rd(pdx, n1) = abs(rd(pdx, n1)) * rat(ro(pdx,1));
    rd(pdx, n2) = abs(rd(pdx, n2)) * rat(ro(pdx,2));
    rd(pdx, n3) = abs(rd(pdx, n3)) * rat(ro(pdx,3));
end

[ARI, sim_sim_agree, CRI, Cout] = cluster_bootstrap(rd, cdx, 500, .75)  ;

[COEFF, SCORE] =  pca(sim_sim_agree);
figure; scatter3(SCORE(:,1), SCORE(:,2),SCORE(:,3), 10, 'k');  rotate3d on
set(gca,'XTicklabel',[], 'YTicklabel', [], 'ZTicklabel', []);