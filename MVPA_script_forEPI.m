db_mode = 0; % debugger skips some of the initial steps and instead loads stuff that usually would be computer from a .mat file
if db_mode == 0;
%mask_fn='/Volumes/Aidas_HDD/MRI_data/S1/Analysis/mask.nii';
mask_fn='/Users/aidas_el_cap/Desktop/ofa_ffa_mask.nii';
%%
roi_path = '/Volumes/Aidas_HDD/MRI_data/S1/Analysis/test3_roi.mat'; load(roi_path);% load to roi object


hold on
disp('preparing for plotting')
for i = 1:4
    hold on
    subplot(4,1,i)
    xlabel(sprintf('Session %d',i))
    drawnow
    P = sprintf('/Volumes/Aidas_HDD/MRI_data/S1/MVPA_analysis/S1/functional/Sess%d/f50hz_wradata.nii',i);
    mY = get_marsy(roi, P, 'mean');  % extract data into marsy data object
y = summary_data(mY);% get summary time course(s)
plot(y)
drawnow
end
disp('done')

% %for plotting
% subplot(2,2,1)
% xlabel('Session 1')
% subplot(2,2,2)
% xlabel('Session 2')
% subplot(2,2,3)
% xlabel('Session 3')
% subplot(2,2,4)
% xlabel('Session 4')

%%
%SPM.xY.P(1:3,:) first three lines
target_EPIs = MakeTRs2;

%filtered_data = '/Volumes/Aidas_HDD/MRI_data/S1/functional/Sess1/f50hz_wrdata.nii';

for ln = 1 : length(target_EPIs)
    hold on
    filename = sprintf('/Volumes/Aidas_HDD/MRI_data/S1/MVPA_analysis/S1/functional/Sess%d/f50hz_wradata.nii',target_EPIs(ln,3));
    single_scan = cosmo_fmri_dataset(filename,'mask', mask_fn,'targets',target_EPIs(ln,2), 'chunks', target_EPIs(ln,3), 'volumes',target_EPIs(ln,1));
    
    subplot(4,1,target_EPIs(ln,3))
plot(target_EPIs(ln,1),mean(single_scan.samples),'r*')
    drawnow
if ln == 1
all_scans=single_scan;
else
all_scans = cosmo_stack({all_scans,single_scan});
disp(['stacking scans ' num2str(length(all_scans.samples(:,1))) ' / ' num2str(length(target_EPIs))])
end
end

%    %% delete expanded scans
%    if delete_expanded_scans == 1
%    for p = 1: 4
%        delete(sprintf('/Volumes/Aidas_HDD/MRI_data/S1/functional/sess%d/f50hz_swrdata_*.nii', p))
%    end
%    disp('deleted expanded scans')
%    end
%     
    
 save('MVPA_stacked_scans')
 disp('saved')
end % end of debug if
%%
if db_mode == 1
cd '/Users/aidas_el_cap/Desktop/MVPA/'
load('MVPA_stacked_scans')
disp('loaded MVPA_stacked_scans.mat')
end

disp('nbrhood')
nbrhood = cosmo_spherical_neighborhood(all_scans, 'radius', 3)
disp('choosing measure')
measure=@cosmo_crossvalidation_measure;  % pick to classify
disp('picking struct')
opt=struct();
opt.classifier=@cosmo_classify_lda;
disp('opt.partitions')
opt.partitions=cosmo_nchoosek_partitioner(all_scans,1);
disp('balancing paritions')
opt.partitions_bal = cosmo_balance_partitions(opt.partitions,all_scans); %added to balance partitions, otherwise cosmo doesnt like it
opt.partitions = opt.partitions_bal 
disp('partitions set and balanced, searchlight is up')

% define spherical neighborhood with radius of 3 voxels

% Run the searchlight with a 3 voxel radius
corr_results=cosmo_searchlight(all_scans,nbrhood,measure,opt);
%
corr_results.samples=corr_results.samples-(1/15);
%% save and exit
output_fn=fullfile(pwd,['outFile.nii']);
output_fn1 = fullfile(pwd,[ datestr(datetime) '_MVPA_outFile.nii'])
try
cosmo_map2fmri(corr_results, output_fn1);
catch
    cosmo_map2fmri(corr_results, output_fn);
end
