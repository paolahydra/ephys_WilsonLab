function [stimset] = bilateralPiezoCommands_windPlayback_ang(iWind,angDownsampling, N_withinPseudoBlock,N_nPseudoBlocks,factor_jointZero,singleStimDur, sample_fs, experimentHandle, armL, armR )
% negative displacements is always pull!
%
% inputs:
%
% iWind has two fields (L, R) and corresponds to one specific wind intensity of
% a fly model (select it and feed it to the function). iWind is sampled
% every 5 degrees (length of each L/R vector: 61)

% angDownsampling, e.g.: 3, to sample every 15 degrees only

% N_withinPseudoBlock = 4;
% factor_jointZero = 2;     % integer. Set to 1 to have the same number of baselines as all the other stimuli. 2 for two as many, and so forth.
% N_nPseudoBlocks = 2;      % with the ability to stop earlier if need so
% singleStimDur = 1.25;     % seconds
% sample_fs = 4e4; 
%

downSamplingVector = 1: angDownsampling: length(iWind.L);
displacementsL = cat(1, zeros(factor_jointZero, 1), armL * tand(iWind.L(downSamplingVector)));
displacementsR = cat(1, zeros(factor_jointZero, 1), armR * tand(iWind.R(downSamplingVector))); %ok

%include so many baselines
% TO BE FINISHED, MAKING SURE I HAVE SOME INDICES TO EASILY SORT THINGS OUT



seconds_addPre = 15;
seconds_addPost = 60;
window_size = 50; % for smoothing
piezo90_ratio_umperVolt = 9;


% voltage_conditions = ang_displacements ./ piezo90_ratio_umperVolt;
N_conditionsTotal = length(downSamplingVector) +1; %adding zero


singleStimDuration_fsp  = sample_fs*singleStimDur;
baseline2add_start_fsp  = sample_fs*seconds_addPre;
baseline2add_end_fsp    = sample_fs*seconds_addPost;


%% fully randomize and build backbone trial-structure

% pre-build
trialTypes = 1 : N_conditionsTotal; %will decode rA and lA from here. No redundancies for now
jointZero_condition = 1; %by design


% build backbone
both = [];
for n = 1: N_nPseudoBlocks
    both_t = repmat([trialTypes, repmat(jointZero_condition, 1, factor_jointZero-1)], 1, N_withinPseudoBlock);
    both = cat(2, both, both_t(randperm(length(both_t))) );
end

duration = length(both)*(singleStimDur); %seconds
duration = duration + seconds_addPre + seconds_addPost;

fprintf('------------\n')
fprintf('n conditions per antenna: %d\n',N_conditionsTotal)
fprintf('n repetitions: %d\n',N_withinPseudoBlock*N_nPseudoBlocks)
fprintf('single sim duration: %.2f s\n',singleStimDur)

fprintf('total duration would be: %d:%02d (min:sec) == %d sec\n', floor(duration/60), round(mod(duration, 60)), round(duration) )
% tabulate(both)
% counts = tabulate_precedingStimulus(both);


nTrials = length(both);
allOnsetsPositions = (1 : singleStimDuration_fsp : nTrials*singleStimDuration_fsp ) + baseline2add_start_fsp; % all trials, independent of type
allOffsetsPositions = (singleStimDuration_fsp : singleStimDuration_fsp : nTrials*singleStimDuration_fsp ) + baseline2add_start_fsp;


%% split backbone into left/right directions
lAnt = both; %legacy
rAnt = both; %legacy

lDisplacements = displacementsL(lAnt);
rDisplacements = displacementsR(rAnt);



%% sample - build command based on DAQ proper. 
% Since it will be background acquisition, I will have to save this to a
% dat file and read and queue just a few rows per time


lCommand = [];
rCommand = [];  % will be inverted, since this probe comes from posterior, so that stim meaning is constant

for n = 1 : length(both)

    lCommand_build = (lDisplacements(n) / piezo90_ratio_umperVolt) * ones(1, sample_fs*singleStimDur);
    lCommand = cat(2, lCommand, lCommand_build);
    
    rCommand_build = (rDisplacements(n) / piezo90_ratio_umperVolt) * ones(1, sample_fs*singleStimDur);
    rCommand = cat(2, rCommand, rCommand_build);
end

rCommand = -rCommand; %IMPORTANT

% add pre- and post- windows before filtering
rCommand = cat(2, zeros(1, baseline2add_start_fsp), rCommand, zeros(1, baseline2add_end_fsp) );
lCommand = cat(2, zeros(1, baseline2add_start_fsp), lCommand, zeros(1, baseline2add_end_fsp) );


clear lCommand_build rCommand_build

%% smooth rather than filtering
rCommand_f = smoothdata(rCommand, 'lowess', window_size);    %comparable to gaussian, better and faster
lCommand_f = smoothdata(lCommand, 'lowess', window_size);    
for i = 1:11
    rCommand_f = smoothdata(rCommand_f, 'lowess', window_size);
    lCommand_f = smoothdata(lCommand_f, 'lowess', window_size);
end  % much better than a single iteration with bigger window: sharper slope, and smoother edges - also faster

commandLR = cat(1, lCommand_f, rCommand_f);

%% write command to file
% rather than outputting the command to workspace, I need to write it to
% file.
identifierkey = datestr(now, 30);
logname = sprintf('p2comm_%s.bin', identifierkey); % will save in current folder wich should be runfolder

fid = fopen(logname, 'w');
fwrite(fid, commandLR, 'single');
fclose(fid);



% % note that I can read it back, specifying 2 rows and as many columns I
% % want my chunk to have. No need to specify the file position, which
% % continues from where the previous read ended, but only if file is never
% % closed and reopened. Example:
% fid = fopen(logname, 'r+');
% A = fread(fid, [2, chuncklengt] ,'single' );
% % when it reaches the end, A is empty.


% 
% %% let's make a simpler example first
% saveData = [1:20; 101:120];
% fid = fopen(logname, 'w');
% fwrite(fid, saveData, 'single');
% fclose(fid);
% 
% 
% %% try to read back chunks
% fid = fopen(logname, 'r+');
% A = fread(fid, [2,5], 'single') %works indeed
% 
% % repeat to see if continues automatically from previous position?
% A = fread(fid, [2,5], 'single') % it does indeed
% 
% % check if this syntax is equivalent
% fseek(fid,0,'cof');
% A = fread(fid, [2,5], 'single') % it is indeed

% fclose(fid)
%% metadata to workspace and saved to disk
stimset.key.ID = identifierkey;
stimset.key.logname = logname;

stimset.joint.trialIndices = both;
stimset.joint.allOnsetsPositions = allOnsetsPositions;
stimset.joint.allOffsetsPositions = allOffsetsPositions;

stimset.singleAnt.L_indices = lAnt;
stimset.singleAnt.R_indices = rAnt;
stimset.singleAnt.L_displacements = lDisplacements; %positive is always push toward head, negative pull
stimset.singleAnt.R_displacements = rDisplacements; %positive is always push toward head, negative pull

%stimset.dec contains redundant information, that may simplify things
%later, possibly
stimset.dec.jointZero_index = jointZero_condition;
stimset.dec.joint2L_decoder_displ = cat(1, 0, iWind.L(downSamplingVector));
stimset.dec.joint2R_decoder_displ = cat(1, 0, iWind.R(downSamplingVector));
stimset.dec.joint2L_decoder_i = trialTypes;
stimset.dec.joint2R_decoder_i = trialTypes;
stimset.dec.i2displ_decoder = iWind; 

stimset.userinput.N_withinPseudoBlock   = N_withinPseudoBlock;
stimset.userinput.factor_jointZero      = factor_jointZero;
stimset.userinput.N_nPseudoBlocks       = N_nPseudoBlocks;
stimset.userinput.singleStimDur         = singleStimDur;
stimset.userinput.DAQ_fs                = sample_fs;

stimset.experimentHandle    = experimentHandle;   % for completeness of info all in the same place, and because this info is not included in the metadata filename


% save to disk
save(sprintf('metadata_%s.mat', identifierkey), '-struct', 'stimset' );


end

