% spmPath='/Volumes/Aidas_HDD/MRI_data/S1/functional/'
% load('/Volumes/Aidas_HDD/MRI_data/S1/Analysis/SPM.mat')
% mask_fn='/Volumes/Aidas_HDD/MRI_data/S1/Analysis/mask.nii';
% %target_betas = [15:29 36:50 57:71]
% %filenames=SPM.xY.P(target_betas,:); % raw epi's get loaded here, filtered,
% %SPM.xY.P(1:3,:) first three lines
% %saved as separate.nii's
% target_EPIs = MakeTRs2;
% %/Volumes/Aidas_HDD/MRI_data/S1/functional/Sess1/swrdata.nii,1 
% 
% %expand the scans
% %%
% expand = 1
% if expand == 1
%     cur_dir = pwd
% disp('Expanding nii scans')
% for p = 1 : 4
%   swr_path = sprintf('/Volumes/Aidas_HDD/MRI_data/S1/functional/sess%d/', p)
%   cd(swr_path)
% expand_nii_scan('swrdata.nii')
% end
% disp('Done expanding')
% cd(cur_dir)
% end
% 
% %%
% 
% base = '/Volumes/Aidas_HDD/MRI_data/S1/functional/'
% for i = 1 : length(target_EPIs)
%     filenames{i,:} = [base 'sess' num2str(target_EPIs(i,3)) '/swrdata_' num2str(target_EPIs(i,1),'%04i'), '.nii'];
% end
% 
% % for i=1:length(target_betas)
% %     filenames(i,:)=[spmPath 'beta_00' num2str(target_betas(i)) '.nii'];
% % end
% %%
% % cc=0;
% % for run=1:4
% %     for task=1:15
% %         cc=cc+1;
% %         single_scan = cosmo_fmri_dataset(filenames(cc,:),'mask', mask_fn,'targets',task, 'chunks', run)
% %         if cc == 1
% %             all_scans=single_scan;
% %         else
% %             all_scans = cosmo_stack({all_scans,single_scan});
% %         end
% %     end
% % end
% %%
% for ln = 1 : length(target_EPIs)
%     single_scan = cosmo_fmri_dataset(filenames{ln},'mask', mask_fn,'targets',target_EPIs(ln,2), 'chunks', target_EPIs(ln,3));
% if ln == 1
% all_scans=single_scan;
% else
% all_scans = cosmo_stack({all_scans,single_scan});
% 
% disp(['stacking scans ' num2str(length(all_scans.samples(:,1))) ' / ' num2str(length(filenames)) ])
% end
% end
% 
%    %% delete expanded scans
%    del_scans = 0;
%    if del_scans == 1
%    for p = 1: 4
%        delete(sprintf('/Volumes/Aidas_HDD/MRI_data/S1/functional/sess%d/swrdata_*.nii', p))
%    end
%    disp('deleted expanded scans')
%    end
%     
%     
%  save('MVPA_stacked_scans')
%  disp('saved')
%%
cd '/Users/aidas_el_cap/Desktop/MVPA/'
load('MVPA_stacked_scans')
disp('loaded MVPA_stacked_scans.mat')
disp('doing nbrhood')
nbrhood = cosmo_spherical_neighborhood(all_scans, 'radius', 2)
disp('choosing measure')
measure=@cosmo_crossvalidation_measure;  % pick to classify
disp('picking struct')
opt=struct();
opt.classifier=@cosmo_classify_lda;
% opt.unbalanced_partitions_ok = true % added because rounded TRs ideally remove this once fixed. 
% unbalanced_partitions_ok = true
% Try partitions=cosmo_balance_partitions(partitions,...).
disp('opt.partitions')
opt.partitions=cosmo_nchoosek_partitioner(all_scans,1);
disp('balancing paritions')
opt.partitions_bal = cosmo_balance_partitions(opt.partitions,all_scans);
opt.partitions = opt.partitions_bal
disp('partitions set and balanced, searchlight is up')

% define spherical neighborhood with radius of 3 voxels

% Run the searchlight with a 3 voxel radius
corr_results=cosmo_searchlight(all_scans,nbrhood,measure,opt);

% 
corr_results.samples=corr_results.samples-(1/15);
% %%
output_fn=fullfile(pwd,['outFile.nii']);
cosmo_map2fmri(corr_results, output_fn);