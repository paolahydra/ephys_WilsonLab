function pipetteNum = runPatchAttempt_continuousAcquisition(ID, pipetteNum)
% % goals:
% measure and record pipette resistance
% record cell-attached trace?
% measure and record resting potential and seal resistance
%
% need to:
% make sure parameters and conversions are right
% re-figure out the exact workflow


HWsettings = sensor_settings_PP; % to set AO channels consistently

daqreset;
s = daq.createSession('ni');
s.Rate = DAQ_fs;
s.IsContinuous = true;



%%
% start acquisition
parametersValue = questdlg_PPStart;


if parametersValue %record trace for paramersn calculation only
    pipetteNum = parametersValue + 1;
    fidWrite = fopen(sprintf('parameters_%s_pipette%02d.bin',ID, pipetteNum), 'w+');
    ai = s.addAnalogInputChannel(HWsettings.DAQdevName, ...
        0:5 , 'Voltage');
    for ch = 1:length(ai)
        ai(ch).Range = [-10 10];
        ai(ch).TerminalConfig = 'SingleEnded';
    end
    lh_ai = s.addlistener('DataAvailable', @(src, event)logData_2Piezos(src, event, fidWrite));
    
    s.startBackground;
    
    stopValue = questdlg_PPStop; %this should work, but I am not sure. Double-check
    s.stop;
    fclose(fidWrite); 
    % now either repeat recursively or record experiment
     
else %record actual experiment
    
    
end

%% record for some time then push a button to stop
% % pause()
% s.stop;
% 
% 
% fclose(fidWrite); 
 





% for subsequent analysis, take from here
%% test data acquis and pip resist
% %read data in
% fid = fopen(sprintf('sealResistanceData_%s.bin',metadata.key.ID), 'r');
% data = fread(fid, 'double' );
% fclose(fid);
% 
% data = reshape(data, length(ai)+1, []);
% timestamps = (data(1,:))';
% data = data(2:end,:);
% data = data';  %for consistency with TO's code
% 
% figure; plot(timestamps(1:4001), data( 1:4001, HWsettings.DATAch.VOLTAGE_NON_SCALED_A_DAQ_AI ))
% hold on
% 
% figure; hold on
% plot(timestamps(1:4001), scaled(1:4001)); ylabel('scaled, [pA]')
% plot(timestamps(1:4001), current(1:4001));
% 
% 
% 
% 
% %% measure pipette resistance
% s.DurationInSeconds = 1;
% [trial_data, trial_time] = s.startForeground();
% [ current, voltage, scaled ] = get_scaled_voltage_and_current_TO( trial_data );
% pipette_resistance = pipetteResistanceCalc( current);
% 
% fprintf('pipette resistance is %3.1f MOhm\n',pipette_resistance)
% 
% 
% %% record access resistance
% % settings.sampRate = 40000; % higher sampling rate for getting the peak
% s.DurationInSeconds = 1.0;
% [trial_data3, trial_time3] = s.startForeground();
% [ current, voltage, scaled ] = get_scaled_voltage_and_current_TO( trial_data3 );
% [R_access, R_membrane] = accessResistanceCalc_TO( current, voltage, trial_time, Ch.sampRate );
% fprintf('access resistance is %3.1f MOhm\n',R_access)
% fprintf('membrane resistance is %3.1f MOhm\n',R_membrane)
% 
% 
% 
% %%  record_i_0_resting_voltage %what was this????
% s.DurationInSeconds = 1.0;
% [trial_data2, trial_time2] = s.startForeground();
% [ current, voltage, scaled ] = get_scaled_voltage_and_current_TO( trial_data2 );
% resting_voltage = mean( voltage );
% fprintf('resting voltage is %3.1f MOhm\n', resting_voltage)
% 
% 
% %%
% daqreset
% 
% save(sprintf('parameters_%s.mat',metadata.key.ID), 'pipette_resistance', 'R_access', 'R_membrane', 'resting_voltage')
