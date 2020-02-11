# HCP_cluster_analysis

All analysis is performed in Matlab (2017a) unless otherwise noted. 

Neuroimaging data were extracted from HCP S900 release, placed in specific folders and organization on my drive. 

Data can be extracted by running:
hcp_clus_initialize.m 
(loads data from tmaps into matlab structure, see file for details)

Data from tasks is stored in the following orders, and extracted from the following contrasts of interest:


Modality:      Con number    Contrast description
1. EMOTION       3           FACES-SHAPES
2. GAMBLING      6           REWARD-PUNISH
3. LANGUAGE      4           STORY-MATH
4. RELATIONAL    4           REL-MATCH
5. SOCIAL        6           TOM-RANDOM
6. WM            11          2BK-0BK

Cognitive data of interest was placed into a CSV file (Extracted from unrestricted_chawco_9_14_2017_11_40_10.csv, not made available due to HCP policies; data can be downloaded from the Human Connectome Project with permission)

The extract cognitive variables, run:
hcp_cog_init.m 

After running these two functions we now have:
hcp_data8.mat (tmaps from HCP data smoothed at 8mm)
hcp_hiearechical.mat (results of hierarchical clustering on all participants)
hcp_cog_dat.mat (cognitive scores and fMRI task performance)

-------

Once cluster membership is calculated, group maps for each cluster was run via:
hcp_group_clus

This function has a dependency on my SPM batch scripts, https://github.com/colinhawco/SPM_bat_scripts

This script ran a group analysis, based on modality (ordered 1-6, see above), with a set of HCP subject IDs, and a clustering solution. This was used in figure 1, and supplemental figures 1-8.

-------

 
cognitive_anovas.m
load clustering data (Wards) and cognitive data do one-way anovas for differences in cognition across clusters, for solutions ranging from 2 to 10
Runs for cognitive test scores, as well as task performance. All results FDR corrected.
Save as cognitive_anovas.mat

This script plots the data in Figure 2 of the paper

-------

The ‘agreement matrices’ (FIGURE 3) were run using the function cluster_bootstrap.m. It runs 1000 permutations of 75% of the data, rerunning the cluster solutions, and generating a probability of any two cases being in the same cluster (from 0 to 1) if both were included in the permutation. It is run from cluster solutions k from 2 to 10. 
Also calculated adjusted Rand index for each cluster solutions (using common participants) for each pair of clusters. This tells us how much overall common information/clustering is present in any two permutations. 
To produce this data I looped through the cluster bootstrap function for modalities and values of k. This is a very computationally demanding process, results are saved to wards_agree.mat

Variables:
Agree: The percentage that any two participants were clustered together, from 0 to 1. The matrix is n x n x k x task, n=822 (sample size), k=cluster number, task =1:6 for the six tasks.
ARI_boots: the adjusted rand index across pairs of permutations; not to save computation it is an upper triangle matrix, lower triangle is zeros. 999 x 1000 x 10 x 6 for the number of permutations (999x1000), k = 2:10, and six tasks. 
RI_boots: the (not adjusted) Rand Index, as ARI_boots
Cout: cluster membership across permutations (0 means not included in sample); can be used to recreate the above three variables. 

Now, we can start really digging into individual figures!! 


------

The script hierarchical_clus_figs.m goes through each figure in turn, plus some visualizations that never made it into the paper. Figures are saved in separate sections and annotated. 

------










