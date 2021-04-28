%% experiment setup
rootfolder = 'D:\Dropbox (HMS)\p2';
experimentHandle = getNewFlyID(rootfolder);             % lets you choose between last run, new run or new fly

% experimentHandle = getNewFlyID(rootfolder, lastFly)   % optional syntax to override last fly specification

disp(experimentHandle)
% output fields:
% 'flyNum'   
% 'runNum'   
% 'basename' 
% 'flyfolder'
% 'runfolder'


%% prepare piezo stims
voltage_conditions = [-2.2, -1.1, 0, 1.1, 2.2];
N_withinPseudoBlock = 4;
factor_jointZero = 2;  % integer. Set to 1 to have the same number of baselines as all the other stimuli. 2 for two as many, and so forth.
N_nPseudoBlocks = 5; %with the ability to stop earlier if need so
singleStimDur = 1.25; %seconds

[stimset,command] = bilateralPiezoCommands(voltage_conditions,N_withinPseudoBlock,N_nPseudoBlocks,factor_jointZero,singleStimDur);



%% setup acquisition
% devices = daq.getDevices;

%settings
Ch.devID = 'Dev1';
Ch.sampRate = 4e4;
ai_channels_used = 0:8;

% set up
s = daq.createSession ('ni');
s.Rate = Ch.sampRate;
aI = s.addAnalogInputChannel(Ch.devID, ai_channels_used, 'Voltage');
for i=1:length(ai_channels_used)
    aI(i).InputType = 'SingleEnded';
end


%% 5 Hz 100mV rectangular waveform - output
% frequency = 1; %Hz
% time_ONorOFF = 0.5/frequency;
% outputVoltage = 0.1; %v to be multiplied by 10, because rear switched input divides by 10
% data = cat(1, ones(time_ONorOFF*Ch.sampRate, 1), zeros(time_ONorOFF*Ch.sampRate, 1)) *outputVoltage *10;
% data = repmat(data, 5, 1);
% 
% aO = s.addAnalogOutputChannel(Ch.devID, 0, 'Voltage');
% s.IsContinuous = true;
% % %
% lh = addlistener(s,'DataRequired', ...
%         @(src,event) src.queueOutputData(data));
% queueOutputData(s,data)
% 
% startBackground(s);
% % %
% % s.stop
% 
% 
% %%
% outputSingleScan(s,1)

%% record pipette resistance (line 400)

s.DurationInSeconds = 1.0;
[trial_data, trial_time] = s.startForeground();
[ current, voltage, scaled ] = get_scaled_voltage_and_current_TO( trial_data );
pipette_resistance = pipetteResistanceCalc( current);

% set(ghandles.pipette_resistance_text, 'String', [num2str(pipette_resistance, 3) ' MOhm']);
fprintf('pipette resistance is %3.1f MOhm\n',pipette_resistance)

% SAVE:

% patch_id = str2num(get(ghandles.patch_id_edit, 'String'));
% experiment_dir = handles.experiment_dir;

% pr_filename = [experiment_dir '/pipette_resistance_patch_A_' num2str(patch_id) '.txt'];
% fileID = fopen(pr_filename,'w+');
% fprintf(fileID,'%f MOhm\n',pipette_resistance);
% fclose(fileID);
% 
% pr_filename_mat = [experiment_dir '/pipette_resistance_patch_A_' num2str(patch_id) '.mat'];
% save(pr_filename_mat, 'current', 'trial_time');



%% record access resistance
% settings.sampRate = 40000; % higher sampling rate for getting the peak
s.DurationInSeconds = 1.0;
[trial_data3, trial_time3] = s.startForeground();
[ current, voltage, scaled ] = get_scaled_voltage_and_current_TO( trial_data3 );
[R_access, R_membrane] = accessResistanceCalc_TO( current, voltage, trial_time, Ch.sampRate );
fprintf('access resistance is %3.1f MOhm\n',R_access)
fprintf('membrane resistance is %3.1f MOhm\n',R_membrane)

%%  record_i_0_resting_voltage
s.DurationInSeconds = 1.0;
[trial_data2, trial_time2] = s.startForeground();
[ current, voltage, scaled ] = get_scaled_voltage_and_current_TO( trial_data2 );
resting_voltage = mean( voltage );
fprintf('resting voltage is %3.1f MOhm\n', resting_voltage)








%% record - change to background acquisition -- or anyways: adapt
s.DurationInSeconds = 60.0;
[trial_dataR, trial_timeR] = s.startForeground();
%  

settings = sensor_settings_PP;
d = regexpdir(experimentHandle.runfolder, 'rec_\w*\d{3}', 0, 0, 0);

save(fullfile(experimentHandle.runfolder, sprintf('rec_%s_%03d.mat', experimentHandle.basename, length(d)+1 )), 'trial_dataR', 'trial_timeR', 'settings' )




