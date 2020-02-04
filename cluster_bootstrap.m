function [ARI, perct_agree, CRI, Cout] = cluster_bootstrap(data, numclus, numboot, bootsize)
% [ARI, perct_agree, CRI, Cout] = cluster_bootstrap(data, numclus, numboot, bootsize)
% bootstrap assessment of cluster solutions
%

if nargin < 4
    bootsize = 0.75;
    if nargin < 3
        numboot = 1000;
    end
end

Z = linkage(data, 'ward');
C = cluster(Z,'MaxClust',numclus);

% establish some variables for later use

nsubs = size(data,1);
nn = 1:nsubs;

% bn is how many subs are in each bootstrap)
bn = ceil(nsubs * bootsize);
perms = zeros(1, nsubs);

agreement=zeros(nsubs);
num_perms = zeros(nsubs);

cout = zeros(nsubs,numboot );

for idx = 1:numboot
    disp(['Creating permuted groupings; Boostrap = ' num2str(idx)])
    test=0;
    % This section makes sure each permutstion of the data is unique.
    while test == 0
        test=1;
        rnd = randperm(nsubs);
        inc(1, 1:nsubs) = 0;
        inc(rnd(1:bn)) = 1;
        % determine if this permutation has been used before
        dist = pdist([perms; inc], 'euclidean');
        if min(dist) == 0
            test=0
        else
            perms(idx,:) = inc;
        end
    end
end

parfor idx = 1:numboot
    disp(['Boostrap = ' num2str(idx)])
    Z = linkage(data(perms(idx,:) == 1, :), 'ward');
    c = cluster(Z,'MaxClust',numclus);
   
    cn=zeros(1,nsubs);
    cn(perms(idx,:)==1) = c;
    cout(:,idx) = cn';
end

Cout = cout;

for idx = 1:numboot-1
    for jdx=idx+1:numboot
        f=find(cout(:,idx) > 0 & cout(:,jdx) > 0);
        c1= cout(f,idx);
        c2=cout(f,jdx);
        [AR,RI,MI,HI]=RandIndex(c1,c2);
        ARI(idx,jdx) = AR;
        CRI(idx,jdx) = RI;
    end
end

% fix the agreement matrix
totagree(1:nsubs,1:nsubs) = 0;  numinc(1:nsubs,1:nsubs) = 0;
for perm = 1:numboot
    dat = cout(:,perm);
    for idx = 1:nsubs
        if dat(idx) > 0;
            match = find(dat == dat(idx));
            totagree(idx,match) = totagree(idx,match)+1;
            inc = find(dat > 0);
            numinc(idx,inc) = numinc(idx,inc)+1;
        end
    end
end
perct_agree =totagree ./ numinc;



