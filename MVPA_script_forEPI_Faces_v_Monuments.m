
%db_mode = 0; % debugger skips some of the initial steps and instead loads stuff that usually would be computer from a .mat file
%if db_mode == 0;

%% Parameters
subID = 5 % which subject
nsess = 5; % how many sessions


plotting = 0
%% Directories
%mask_fn='/Volumes/Aidas_HDD/MRI_data/S1/Analysis/mask.nii';
%mask_fn= sprintf('/Volumes/Aidas_HDD/MRI_data/S%d/Analysis/mask.nii',subID); % brain mask used for searchlight 
mask_fn = '/Volumes/Aidas_HDD/MRI_data/MVPA_analyses/MVP_mask_for_MVPA.nii' % ROI-ish mask
addpath('/Users/aidas_el_cap/Documents/MATLAB/spm12/toolbox/marsbar/'); %marsbar dir keeps disappearing re-add it just in case
mvpadir = '/Volumes/Aidas_HDD/MRI_data/MVPA_analyses/'; % where to save stacked scans and outpt files
%% Load myTrials, plot and prep
load(sprintf('/Volumes/Aidas_HDD/MRI_data/Other/myTrials3456_processed/S%d_Results.mat',subID));
myTrials = MakeTRs2_facesvmonuments(myTrials);
%% faces = target #1, monuments = target #2
ffc = find([myTrials.name_ID] > 0);
[myTrials(ffc).name_ID] = deal(1);
%%
% Done preparing and plotting at this point.

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
%target_EPIs =;
%load('/Volumes/Aidas_HDD/MRI_data/S2/fixed_target_EPIs.mat')

%filtered_data = '/Volumes/Aidas_HDD/MRI_data/S1/functional/Sess1/f50hz_wrdata.nii';
%scans = find([myTrials.name_ID] ~= 0);
%scans = find([myTrials.blockNum] ~= 14);
scans = find([myTrials.blockNum] == 15 | [myTrials.blockNum] == 16 | [myTrials.blockNum] == 17 | [myTrials.blockNum] == 1 | [myTrials.blockNum] == 2 | [myTrials.blockNum] == 3)
%scans = 1 : length(myTrials);
%%
for ln = scans
    disp(ln)
    %hold on
    %filename = sprintf('/Volumes/Aidas_HDD/MRI_data/S2/Functional/Sess%d/f50hz_swradata.nii',target_EPIs(ln,3));
    filename = sprintf('/Volumes/Aidas_HDD/MRI_data/S%d/Functional/Sess%d/f50hz_swradata.nii',subID,myTrials(ln).fmriRun);
    single_scan = cosmo_fmri_dataset(filename,'mask', mask_fn,'targets',myTrials(ln).name_ID, 'chunks', myTrials(ln).fmriRun, 'volumes',myTrials(ln).TR);
    %subplot(4,1,myTrials(ln).fmriRun)
    %plot(myTrials(ln).TR,mean(single_scan.samples),'r*')
    %drawnow
  if ln == 1
    all_scans=single_scan;
    else
    all_scans = cosmo_stack({all_scans,single_scan});
  end
  single_scan = cosmo_fmri_dataset(filename,'mask', mask_fn,'targets',myTrials(ln).name_ID, 'chunks', myTrials(ln).fmriRun, 'volumes', myTrials(ln).TR + 1);
  all_scans = cosmo_stack({all_scans,single_scan});
  %ends ln == 1 if statement
    %disp(['stacking scans ' num2str(length(all_scans.samples(:,1))) ' / ' num2str(length(scans)
end
%%
disp('saving stacked scans')
save(fullfile(mvpadir,sprintf('MVPA_stacked_scansVmon_sub%d.mat',subID)),'all_scans')
disp('saved')

%% Z Scoring
% 1 way, grab task in every run, zscore
raw_all_scans = all_scans; %
tasks = unique(all_scans.sa.chunks);
for p = 1 : length(tasks);
tinx = find(all_scans.sa.chunks == p);
all_scans.samples(tinx,:) = zscore(all_scans.samples(tinx,:),[],1);
end

%%
%    %% delete expanded scans
%    if delete_expanded_scans == 1
%    for p = 1: 4
%        delete(sprintf('/Volumes/Aidas_HDD/MRI_data/S1/functional/sess%d/f50hz_swrdata_*.nii', p))
%    end
%    disp('deleted expanded scans')
%    end
%date_nao = datestr(date); lol weirdest error

%exist(sprintf('MVPA_stacked_scans_sub%d.mat',subID),'file') %checks if
%stacked scans already exist
 
%  end%end of debuf if
%%
% if db_mode == 1
% cd '/Users/aidas_el_cap/Desktop/MVPA/'
% load('MVPA_stacked_scans')
% disp('loaded MVPA_stacked_scans.mat')
% end
%%
disp('nbrhood')
nbrhood = cosmo_spherical_neighborhood(all_scans, 'radius', 3)
disp('choosing measure')
measure=@cosmo_crossvalidation_measure;  % pick to classify
disp('picking struct')
opt=struct();
opt.classifier=@cosmo_classify_lda;
disp('opt.partitions')
%opt.partitions=cosmo_nchoosek_partitioner(all_scans,1);
opt.partitions = cosmo_nfold_partitioner(all_scans);
%% Balance
disp('balancing paritions')
opt.partitions_bal = cosmo_balance_partitions(opt.partitions,all_scans,'nrepeats',3); %added to balance partitions, otherwise cosmo doesnt like it
opt.partitions = opt.partitions_bal 
disp('partitions set and balanced, searchlight is up')
%%


cosmo_check_dataset(all_scans)
cosmo_check_partitions(opt.partitions,all_scans)

l = find(all_scans.sa.targets == 0)
all_scans.sa.targets(l) = 2
% define spherical neighborhood with radius of 3 voxels

% Run the searchlight with a 3 voxel radius
corr_results=cosmo_searchlight(all_scans,nbrhood,measure,opt); % ERRORS HERE
corr_results.samples=corr_results.samples-(1/40);
%% save and exit
output_fn=fullfile(mvpadir,sprintf('Sub%d_MVPA_results.nii',subID));
cosmo_map2fmri(corr_results, output_fn);
save(sprintf('Sub%d_MVPA_wrkspc',subID))
% output_fn1 = fullfile(pwd,[ datestr(datetime) '_MVPA_outFile.nii'])
% try
% cosmo_map2fmri(corr_results, output_fn1);
% catch
%     cosmo_map2fmri(corr_results, output_fn);

cosmo_check_partitions(opt.partitions,all_scans)
opt.partitions = cosmo_nfold_partitioner(all_scans)
cosmo_balance_partitions(opt.partitions,all_scans)