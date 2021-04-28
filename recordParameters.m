function [pipetteNum, parametersValue] = recordParameters(s, ID, pipetteNum )

    pipetteNum = pipetteNum+1;
    fidWrite = fopen(sprintf('parameters_%s_pipette%02d.bin',ID, pipetteNum), 'w+');
    lh_ai = s.addlistener('DataAvailable', @(src, event)logData_2Piezos(src, event, fidWrite));
    
    s.startBackground;
    parametersValueLocal = questdlg_PPContinue_orStop; % 1 continue, 2 stop
       
    %% calculate pipette resistance, then continue acquisition until stop
    if parametersValueLocal == 1
        s.stop;
        
        % read in the last x seconds
        timeRead = 1; %seconds
        nBits = 8;
        fseek(fidWrite, -(length(s.Channels)+1) * timeRead * s.Rate * nBits, 'eof');
        data = fread(fidWrite, 'double');
        
        data = reshape(data, length(s.Channels)+1, []);
        endTime = data(1,end);
        data = data(2:end,:);
        data = data';  %for consistency with TO's code
        
        [ current, ~, ~ ] = get_scaled_voltage_and_current_PP( data );
        parameters.pipette_resistance = pipetteResistanceCalc( current);
        parameters.pip_resist_endTime = endTime;
        
        fprintf('pipette resistance is %3.1f MOhm\n',parameters.pipette_resistance)
        save(sprintf('parameters_%s_pipette%02d.mat',ID, pipetteNum), 'parameters');
        
        
        % and when you are done resume recording.
        disp('resuming recording...')
        % I don't think I have to add the listener back?
        s.startBackground;
        questdlg_PPStop;
        
    end
    

    s.stop;   
    fclose(fidWrite);
    delete(lh_ai)
    
    parametersValue = questdlg_PPStart; 
    if parametersValue==1 %record trace for paramersn calculation only
        [pipetteNum, parametersValue] = recordParameters(s, ID, pipetteNum ); % pipetteNum updated in here
    end

end

