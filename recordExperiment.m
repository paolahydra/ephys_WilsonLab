function recordExperiment(s, HWsettings, metadata, DAQ_fs, pipetteNum )
    % add input channels
    fidWrite = fopen(sprintf('ephysData_%s_pipette%02d.bin',metadata.key.ID, pipetteNum), 'w+');
    s.addAnalogInputChannel(HWsettings.DAQdevName, [HWsettings.AIch.PIEZO_SENSOR_LEFT, HWsettings.AIch.DUMMY1, HWsettings.AIch.PIEZO_SENSOR_RIGHT] , 'Voltage');
    for ch = 7:9
        s.Channels(ch).Range = [-10 10];
        s.Channels(ch).TerminalConfig = 'SingleEnded';
%         disp('SGS channels set to differential')
    end
    lh_ai = s.addlistener('DataAvailable', @(src, event)logData_2Piezos(src, event, fidWrite)); 
    
    
    
    % add output channels
    fidRead = fopen(metadata.key.logname, 'r+'); %only close it at the end
    chunklength = 12*DAQ_fs;
    s.addAnalogOutputChannel(HWsettings.DAQdevName, ...
        [HWsettings.AOch.PIEZO_OUTPUT_LEFT, HWsettings.AOch.PIEZO_OUTPUT_RIGHT] , 'Voltage');
    for ch = 10:11
        s.Channels(ch).Range = [-10 10];
        s.Channels(ch).TerminalConfig = 'SingleEnded';
    end
    % do the first fread here, and then continue
    outData = fread(fidRead, [2, chunklength] ,'single' );
    queueOutputData(s, outData'); %these need to be one column per ao channel
    lh_ao = s.addlistener('DataRequired', @(src, event)queueMoreData_2Piezos(src, event, fidRead, chunklength));

    
    
    
    %%
    s.startBackground;
    
    tic
    FS = stoploop({'Stop'}) ;
    % Display elapsed time
    fprintf('\nRECORDING: elapsed time (s): %8.2f\t',toc)
    % start the loop
    while(~FS.Stop() && s.IsRunning)       % Check if the loop has to be stopped
        fprintf('%c',repmat(8,9,1)) ;   % clear up previous time
        fprintf('%8.2f\t',toc) ;        % display elapsed time
    end
    
    
    %either finished, or user stopped acquisition
    if FS.Stop() % user stopped
        s.stop;
    end
    FS.Clear() ;  % Clear up the box
    clear FS ;
    
    
    outputSingleScan(s,[0,0]);
    
    
    fclose(fidWrite);
    fclose(fidRead);
    delete(lh_ai);
    delete(lh_ao);  
    delete(s);
    save(sprintf('acqSettings_%s.mat',metadata.key.ID), 'HWsettings')

end

