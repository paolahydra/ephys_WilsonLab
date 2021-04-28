% record either parameters (recursive) or experiment (final), or return


%% general settings
HWsettings = sensor_settings_PP; % to set AO channels consistently
daqreset;
s = daq.createSession('ni');
s.Rate = DAQ_fs;
s.IsContinuous = true;
% s.NotifyWhenDataAvailableExceeds = 2000;    % 4000
% s.NotifyWhenScansQueuedBelow = 20000;       % 20000
s.addAnalogInputChannel(HWsettings.DAQdevName, 0:5 , 'Voltage'); %common across
for ch = 1:length(s.Channels)
    s.Channels(ch).Range = [-10 10];
    s.Channels(ch).TerminalConfig = 'SingleEnded';
end

%%
parametersValue = questdlg_PPStart;     % first time asking

if parametersValue==1 %record trace for paramersn calculation only
   [pipetteNum, parametersValue] = recordParameters(s, metadata.key.ID, pipetteNum ); % pipetteNum updated in here
end

if parametersValue==2 %record actual experiment
   recordExperiment(s, HWsettings, metadata, DAQ_fs, pipetteNum )
end



