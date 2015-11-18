spmPath='/Volumes/Aidas_HDD/MRI_data/S1/Analysis/'
load('/Volumes/Aidas_HDD/MRI_data/S1/Analysis/SPM.mat')
mask_fn='/Volumes/Aidas_HDD/MRI_data/S1/Analysis/mask.nii';
target_betas = [15:29 36:50 57:71]
%filenames=SPM.xY.P(target_betas,:); % raw epi's get loaded here, filtered,
%saved as separate.nii's

for i=1:length(target_betas)
    filenames(i,:)=[spmPath 'beta_00' num2str(target_betas(i)) '.nii'];
end
%%
cc=0;
for run=1:3
    for task=1:15
        cc=cc+1;
        single_scan = cosmo_fmri_dataset(filenames(cc,:),'mask', mask_fn,'targets',task, 'chunks', run)
        if cc == 1
            all_scans=single_scan;
        else
            all_scans = cosmo_stack({all_scans,single_scan});
        end
    end
end
%%

nbrhood = cosmo_spherical_neighborhood(all_scans, 'radius', 3)
measure=@cosmo_crossvalidation_measure;  % pick to classify
opt=struct();
opt.classifier=@cosmo_classify_lda;
opt.partitions=cosmo_nchoosek_partitioner(all_scans,1);

% define spherical neighborhood with radius of 3 voxels

% Run the searchlight with a 3 voxel radius
corr_results=cosmo_searchlight(all_scans,nbrhood,measure,opt);

% 
corr_results.samples=corr_results.samples-(1/15);
% %%
output_fn=fullfile(pwd,['outFile.nii']);
cosmo_map2fmri(corr_results, output_fn);