function [pipetteNum, parametersValue] = recordParameters(s, ID, pipetteNum )

    pipetteNum = pipetteNum+1;
    fidWrite = fopen(sprintf('parameters_%s_pipette%02d.bin',ID, pipetteNum), 'w+');
    lh_ai = s.addlistener('DataAvailable', @(src, event)logData_2Piezos(src, event, fidWrite));
    
    s.startBackground;
    
    questdlg_PPStop; %this should work, but I am not sure. Double-check
    %%
    s.stop;
    fclose(fidWrite);
    delete(lh_ai)
    
    parametersValue = questdlg_PPStart; 
    if parametersValue==1 %record trace for paramersn calculation only
        [pipetteNum, parametersValue] = recordParameters(s, ID, pipetteNum ); % pipetteNum updated in here
    end

end

