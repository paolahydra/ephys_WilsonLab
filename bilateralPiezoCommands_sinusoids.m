function [stimset] = bilateralPiezoCommands_sinusoids(displacements, interleaveVibrations, N_withinPseudoBlock,N_nPseudoBlocks,factor_jointZero,singleStimDur, sample_fs, experimentHandle )
% negative displacements is always pull!
%
% inputs:
%
N_withinPseudoBlock = 1;
factor_jointZero = 2;     % integer. Set to 1 to have the same number of baselines as all the other stimuli. 2 for two as many, and so forth.
N_nPseudoBlocks = 20;      % with the ability to stop earlier if need so
singleStimDur = 1;     % seconds
sample_fs = 4e4; 
window_size = 50; % for smoothing
%
seconds_addPre = 15;
seconds_addPost = 15;
piezo90_ratio_umperVolt = 9;


displacements = [-28   -14    0    14   28];
interleaveVibrations = 1;



% hard-build vibration carrierFr and Amplitude based on displacement levels
% used
for d = 1:length(displacements)
    switch abs(displacements(d))
        case 28
            VibrCarrierFr(d)    = 400;
            VibrAmplitudeUM(d)  = 2; %need piezo-specific correction for freq dependent amplitude drop
        case 14
            VibrCarrierFr(d)    = 250;
            VibrAmplitudeUM(d)  = 1; %need piezo-specific correction for freq dependent amplitude drop
%         case 7
%             VibrCarrierFr(d)    = 275;
%             vibrAmplitudeUM(d)  = 0.5; %need piezo-specific correction for freq dependent amplitude drop
        case 0
            VibrCarrierFr(d)    = 100;
            VibrAmplitudeUM(d)  = 2.5; %need piezo-specific correction for freq dependent amplitude drop
    end
end
VibrAmplitudeVolts = VibrAmplitudeUM ./ piezo90_ratio_umperVolt;
correctFactorL = [1.3,  1.1,    1,      1.1,    1.3  ]; %ones(size(VibrAmplitudeVolts));
correctFactorR = [7,    3.2,  1.4,      3.2,    7  ]; %ones(size(VibrAmplitudeVolts));


voltage_conditions = displacements ./ piezo90_ratio_umperVolt;
N_conditionsPerAntenna = length(voltage_conditions);


singleStimDuration_fsp  = sample_fs*singleStimDur;
baseline2add_start_fsp  = sample_fs*seconds_addPre;
baseline2add_end_fsp    = sample_fs*seconds_addPost;


%object with shared properties
vib = PipStimulus;
vib.modulationFreqHz = 2;
vib.modulationPhase = 0.25; % -0.25 starts low, 0.25 starts high
vib.pipDur = singleStimDur;
vib.endPadDur = 0;
vib.startPadDur = 0;
vib.sampleRate = sample_fs;
vib.amplitude = 0.1;
vib.maxVoltage = 1;

%% fully randomize and build backbone trial-structure

% pre-build
trialTypes = 1 : N_conditionsPerAntenna^2; %will decode rA and lA from here. No redundancies for now
[L,R] = ind2sub([N_conditionsPerAntenna, N_conditionsPerAntenna],trialTypes); %this is used to map indices. No redundancies in here
% to overrepresent joint baseline:
i_vZero = find(voltage_conditions==0);
jointZero_condition = find(L==i_vZero & R==i_vZero);


% build backbone
both = [];      % now by columns
for n = 1: N_nPseudoBlocks
    both_t = repmat([trialTypes, repmat(jointZero_condition, 1, factor_jointZero-1)], 1, N_withinPseudoBlock);
    both_t = both_t';
    both = cat(2, both, both_t(randperm(length(both_t))) );
end
clear both_t

if interleaveVibrations
    both = cat(1, both, both+N_conditionsPerAntenna^2);  % trialTypes>N_conditionsPerAntenna^2 represent vibration trials
end
both = both(:)';



duration = length(both)*(singleStimDur); %seconds
duration = duration + seconds_addPre + seconds_addPost;

fprintf('------------\n')
if interleaveVibrations
    fprintf('n conditions per antenna: %d\n',N_conditionsPerAntenna*2)
else
    fprintf('n conditions per antenna: %d\n',N_conditionsPerAntenna)
end
fprintf('n repetitions: %d\n',N_withinPseudoBlock*N_nPseudoBlocks)
fprintf('single sim duration: %.2f s\n',singleStimDur)

fprintf('total duration would be: %d:%02d (min:sec) == %d sec\n', floor(duration/60), round(mod(duration, 60)), round(duration) )
% tabulate(both)
% counts = tabulate_precedingStimulus(both);


nTrials = length(both);
allOnsetsPositions = (1 : singleStimDuration_fsp : nTrials*singleStimDuration_fsp ) + baseline2add_start_fsp; % all trials, independent of type
allOffsetsPositions = (singleStimDuration_fsp : singleStimDuration_fsp : nTrials*singleStimDuration_fsp ) + baseline2add_start_fsp;


%% split backbone into left/right directions - sustained only
for n = 1 : length(both)
    sustStim = mod(both(n), N_conditionsPerAntenna^2 );
    if sustStim == 0
        sustStim =  N_conditionsPerAntenna^2 ;
    end
    lAnt(n) = L(sustStim);
    rAnt(n) = R(sustStim);
    
    lDisplacements(n) = displacements(lAnt(n));
    rDisplacements(n) = displacements(rAnt(n));
end
clear sustStim

%% sample - build command based on DAQ proper. 
% Since it will be background acquisition, I will have to save this to a
% dat file and read and queue just a few rows per time

lCommand = [];
rCommand = [];  % will be inverted, since this probe comes from posterior, so that stim meaning is constant
lCommandVibr = [];
rCommandVibr = [];

for n = 1 : length(both)
    disp(n)
    lCommand_build = voltage_conditions(lAnt(n)) * ones(1, sample_fs*singleStimDur);
    lCommandVibr_build = zeros(size(lCommand_build));
    if both(n) > N_conditionsPerAntenna^2
        % produce stimulus-specific vibration
        vib.carrierFreqHz = VibrCarrierFr(lAnt(n));
        vib.amplitude = VibrAmplitudeVolts(lAnt(n)) * correctFactorL(lAnt(n)) /10;
        lCommandVibr_build = (vib.stimulus)';
    end
    lCommand = cat(2, lCommand, lCommand_build);
    lCommandVibr = cat(2, lCommandVibr, lCommandVibr_build);
    
    
    rCommand_build = voltage_conditions(rAnt(n)) * ones(1, sample_fs*singleStimDur);
    rCommandVibr_build = zeros(size(rCommand_build));
    if both(n) > N_conditionsPerAntenna^2
        %add stimulus-specific vibration
        vib.carrierFreqHz = VibrCarrierFr(rAnt(n));
        vib.amplitude = VibrAmplitudeVolts(rAnt(n)) * correctFactorR(rAnt(n)) /10;
        rCommandVibr_build = (vib.stimulus)';
    end
    rCommand = cat(2, rCommand, rCommand_build);
    rCommandVibr = cat(2, rCommandVibr, rCommandVibr_build);
end

rCommand = -rCommand; %IMPORTANT

% add pre- and post- windows before filtering
rCommand = cat(2, zeros(1, baseline2add_start_fsp), rCommand, zeros(1, baseline2add_end_fsp) );
lCommand = cat(2, zeros(1, baseline2add_start_fsp), lCommand, zeros(1, baseline2add_end_fsp) );

lCommandVibr = cat(2, zeros(1, baseline2add_start_fsp), lCommandVibr, zeros(1, baseline2add_end_fsp) );
rCommandVibr = cat(2, zeros(1, baseline2add_start_fsp), rCommandVibr, zeros(1, baseline2add_end_fsp) );

lCommandVibr =lCommandVibr.* 10;
rCommandVibr =rCommandVibr.* 10;

clear lCommand_build rCommand_build

%% smooth rather than filtering
rCommand_f = smoothdata(rCommand, 'lowess', window_size);    %comparable to gaussian, better and faster
lCommand_f = smoothdata(lCommand, 'lowess', window_size);    
for i = 1:11
    rCommand_f = smoothdata(rCommand_f, 'lowess', window_size);
    lCommand_f = smoothdata(lCommand_f, 'lowess', window_size);
end  % much better than a single iteration with bigger window: sharper slope, and smoother edges - also faster


commandLR = cat(1, lCommand_f+lCommandVibr, rCommand_f+rCommandVibr);
% figure; plot(commandLR(1,:)); hold on; plot(commandLR(2,:))
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
stimset.dec.joint2L_decoder_displ = displacements(L);
stimset.dec.joint2R_decoder_displ = displacements(R);
stimset.dec.joint2L_decoder_i = L;
stimset.dec.joint2R_decoder_i = R;
stimset.dec.i2displ_decoder = displacements;

stimset.userinput.N_withinPseudoBlock   = N_withinPseudoBlock;
stimset.userinput.factor_jointZero      = factor_jointZero;
stimset.userinput.N_nPseudoBlocks       = N_nPseudoBlocks;
stimset.userinput.singleStimDur         = singleStimDur;
stimset.userinput.DAQ_fs                = sample_fs;

stimset.experimentHandle    = experimentHandle;   % for completeness of info all in the same place, and because this info is not included in the metadata filename


% save to disk
save(sprintf('metadata_%s.mat', identifierkey), '-struct', 'stimset' );


end

