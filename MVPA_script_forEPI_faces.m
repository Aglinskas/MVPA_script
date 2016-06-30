description = 'real per block z scoring only good faces, only single scan loaded, 6 second hrf offset'
%% Parameters 
subID = 5 % which subject
nsess = 5; % how many sessions
plotting = 0 % enable plotting?
%% Directories % filenames
%mask_fn='/Volumes/Aidas_HDD/MRI_data/S1/Analysis/mask.nii';
%mask_fn= sprintf('/Volumes/Aidas_HDD/MRI_data/S%d/Analysis/mask.nii',subID); % brain mask used for searchlight 
%mask_fn = '/Volumes/Aidas_HDD/MRI_data/MVPA_analyses/MVP_mask_for_MVPA.nii' % ROI-ish mask
%mask_fn = '/Volumes/Aidas_HDD/MRI_data/S4/Analysis2/mask.nii'
mask_fn = sprintf('/Volumes/Aidas_HDD/MRI_data/S%d/Analysis2/mask.nii',subID) %Brain mask from GLM
%mask_fn = '/Users/aidas_el_cap/Desktop/MEGA_ROI_s4.nii'
%addpath('/Users/aidas_el_cap/Documents/MATLAB/spm12/toolbox/marsbar/'); %marsbar dir keeps disappearing re-add it just in case
mvpadir = '/Volumes/Aidas_HDD/MRI_data/MVPA_analyses/'; % where to save stacked scans and outpt files

%% Load myTrials, plot and prep
load(sprintf('/Volumes/Aidas_HDD/MRI_data/Other/myTrials3456_processed/S%d_Results.mat',subID));
myTrials = MakeTRs2_faces(myTrials);
% if exist(fullfile(mvpadir,sprintf('MVPA_stacked_scans_sub%d.mat',subID)),'file') == 2
%     disp('Stacked scans found, skipping plotting and stacking')
% end
% if exist(fullfile(mvpadir,sprintf('MVPA_stacked_scans_sub%d.mat',subID)),'file') ~= 2
% 
% 
% 
% try
% roi_path = sprintf('/Volumes/Aidas_HDD/MRI_data/S%d/Analysis/FFA_ROI_roi.mat',subID); load(roi_path); %FFAish_roi.mat'; load(roi_path);% load to roi object
% end
% 
% if plotting == 1
% clf;
% hold on;
% disp('preparing for plotting')
% for i = 1:nsess
%     %hold on
%     subplot(nsess,1,i)
%     xlabel(sprintf('Session %d',i))
%     drawnow
%     P = sprintf('/Volumes/Aidas_HDD/MRI_data/S%d/Functional/Sess%d/f50hz_rrdata.nii',subID,i);
%     mY = get_marsy(roi, P, 'mean'); % extract data into marsy data object   
%     % get_marsy_1 <- pure magic 
% y = summary_data(mY);% get summary time course(s) % needs marsy object as input
% plot(y)
% drawnow
% end
% disp('done')
% end
% end
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
% if exist(fullfile(mvpadir,sprintf('MVPA_stacked_scans_sub%d.mat',subID)),'file') == 2 % if stacked scans exist
%     load(fullfile(mvpadir,sprintf('MVPA_stacked_scans_sub%d.mat',subID)))
%     disp('previously stacked scans loaded')
% else 
%     disp('stacking scans')
%scans = find([myTrials.name_ID] ~= 0);
%%
% fix skipped answers
for i = 1 : length(myTrials)
if isempty(myTrials(i).resp)
myTrials(i).resp = 0
end
end
%% 
% bad faces are: 
% if they answered no to "7. do you know this person", or 
%'I don't don't know this person' to question 11. "what does this person do?"'

bad_faces = [[myTrials(find([myTrials.blockNum] == 7 & [myTrials.resp] == 2)).name_ID]',
[myTrials(find([myTrials.blockNum] == 11 & [myTrials.resp] == 4)).name_ID]']
% [myTrials(find([myTrials.blockNum] == 14 & [myTrials.resp] == 4)).name_ID]'];
% [myTrials(find([myTrials.blockNum] == 2 & [myTrials.resp] == 4)).name_ID]'
% [myTrials(find([myTrials.blockNum] == 9 & [myTrials.resp] == 4)).name_ID]'];


%legacy
%scans = find([myTrials.name_ID] ~= 0 & [myTrials.blockNum] ~= 14);
%only good faces
scans = find([myTrials.name_ID] ~= 0 & [myTrials.blockNum] ~= 14 & ismember([myTrials.name_ID],bad_faces) ~= 1);
% %%
for ln = scans
    disp(ln)
    %hold on
    %filename = sprintf('/Volumes/Aidas_HDD/MRI_data/S2/Functional/Sess%d/f50hz_swradata.nii',target_EPIs(ln,3));
    %filename = sprintf('/Volumes/Aidas_HDD/MRI_data/S%d/Functional/Sess%d/f50hz_swradata.nii',subID,myTrials(ln).fmriRun);
    filename = sprintf('/Volumes/Aidas_HDD/MRI_data/S%d/Functional/Sess%d/f50hz_rrdata.nii',subID,myTrials(ln).fmriRun);
    single_scan = cosmo_fmri_dataset(filename,'mask', mask_fn,'targets',myTrials(ln).name_ID, 'chunks', myTrials(ln).blockNum, 'volumes',myTrials(ln).TR);
    %subplot(4,1,myTrials(ln).fmriRun)
    %plot(myTrials(ln).TR,mean(single_scan.samples),'r*')
    %drawnow
  if ln == 1
    all_scans=single_scan;
    else
    all_scans = cosmo_stack({all_scans,single_scan});
  end
  %single_scan = cosmo_fmri_dataset(filename,'mask', mask_fn,'targets',myTrials(ln).name_ID, 'chunks', myTrials(ln).blockNum, 'volumes', myTrials(ln).TR + 1);
  %all_scans = cosmo_stack({all_scans,single_scan});
  %ends ln == 1 if statement
    %disp(['stacking scans ' num2str(length(all_scans.samples(:,1))) ' / ' num2str(length(scans)
end
%%
disp('saving stacked scans')
save(fullfile(mvpadir,sprintf('MVPA_stacked_scans_sub%d.mat',subID)),'all_scans')
disp('saved')

% Z Scoring
%first way, grab task in every run, zscore

% raw_all_scans = all_scans; %back of of raw scores;
% tasks = unique(all_scans.sa.chunks); %loads the tasks that are used in chunks
% for p = 1 : length(tasks);
% tinx = find(all_scans.sa.chunks == p); %finds trials belongng to that task in all runs
% all_scans.samples(tinx,:) = zscore(all_scans.samples(tinx,:),[],1);
% end

% %% Same but softmax
% 
% raw_all_scans = all_scans; %back of of raw scores;
% tasks = unique(all_scans.sa.chunks); %loads the tasks that are used in chunks
% for p = 1 : length(tasks);
% tinx = find(all_scans.sa.chunks == p); %finds trials belongng to that task in all runs
% all_scans.samples(tinx,:) = softmax(all_scans.samples(tinx,:),[],1);
% end

%% Per block z scoring, warning: ugly code, will clean up in near future

raw_all_scans = all_scans; % make a backup
i = 0; % init counters
bl_s = 0;
while i < length(all_scans.sa.chunks) %#468 go through all chunks


    bl_s = i + 1;
%bl_s = all_scans.sa.chunks(14);

i = bl_s;


while all_scans.sa.chunks(i + 1) == all_scans.sa.chunks(i); %count the number of scans in the block
   
    %c_bl_length = c_bl_length + 1
    i = i + 1;
%     if i == length(all_scans.sa.chunks) - 1;
%     break
% end 
end

if i == length(all_scans.sa.chunks) - 1; % manualy break out at the end
    break
end 

i;
bl_length = i - bl_s; %block length, not super necessary, but nice to have for post-mortem de-bugging

[bl_s bl_length i];
all_scans.sa.chunks(bl_s,2) = 1;
all_scans.sa.chunks(i,2) = 2;

all_scans.samples(bl_s:i,:) = zscore(all_scans.samples(bl_s:i,:),[],1) % Actual z scoring
end
debug_zScoring = all_scans.sa.chunks; %make a copy of block index
all_scans.sa.chunks(:,2) = [] %
%%
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
opt.partitions=cosmo_nchoosek_partitioner(all_scans,1);

%% other available options:
%opt=struct();
%     opt.output='accuracy';
%     opt.normalization='zscore';
%     opt.classifier=@cosmo_classify_lda;
%     opt.partitions=cosmo_nchoosek_partitioner(ds,2);


%disp('balancing paritions')
%opt.partitions_bal = cosmo_balance_partitions(opt.partitions,all_scans); %added to balance partitions, otherwise cosmo doesnt like it
%opt.partitions = opt.partitions_bal 
%disp('partitions set and balanced, searchlight is up')

% define spherical neighborhood with radius of 3 voxels

% Run the searchlight with a 3 voxel radius
%%
corr_results=cosmo_searchlight(all_scans,nbrhood,measure,opt); %
%corr_results.samples=corr_results.samples-(1/40);
corr_results.samples=corr_results.samples-(1/length(unique(all_scans.sa.targets)));

% r structure has the accuracies
r{1,1} = 'mean accuracy';
r{1,2} = 'max accuracy';
r{1,3} = '#vx above 2% accuracy';
r{1,4} = '#vx above 3% accuracy';
r{1,5} = '#vx above 4% accuracy';
r{2,1} = mean(corr_results.samples);
r{2,2} = max(corr_results.samples);
r{2,3} = length(find(corr_results.samples > 0.02));
r{2,4} = length(find(corr_results.samples > 0.03));
r{2,5} = length(find(corr_results.samples > 0.04));

%% save and exit
file_name = sprintf('Sub%d_MVPA_results',subID)
output_fn=fullfile(mvpadir,file_name);

wh = 1
while exist([output_fn num2str(wh) '.nii'],'file') > 0
    wh = wh + 1
end
cosmo_map2fmri(corr_results,[output_fn num2str(wh) '.nii']);

file_name2 = sprintf('Sub%d_MVPA_wrkspc',subID)
output_fn2 = fullfile(mvpadir,file_name2)
wh1 = 1
while exist([output_fn2 num2str(wh1) '.mat'],'file') > 0
    wh1 = wh1 + 1
end
save([output_fn2 num2str(wh1) '.mat'])



%end
% output_fn1 = fullfile(pwd,[ datestr(datetime) '_MVPA_outFile.nii'])
% try
% cosmo_map2fmri(corr_results, output_fn1);
% catch
%     cosmo_map2fmri(corr_results, output_fn);
