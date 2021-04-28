function parametersValue = questdlg_PPContinue_orStop
    d = dialog('Position',[2951, 1075, 210, 60],'Name','Stop');

    pipResist  = uicontrol('Parent',d,...
               'Position',[17 10 95 40],...
               'String',' PIP RESIST',...
               'Callback',@pipetteResistanceAndContinue_callback);
           
    experimBtn = uicontrol('Parent',d,...
               'Position',[126 10 70 40],...
               'String','STOP',...
               'Callback',@stopParameters_callback);
    
    parametersValue = 0;       
    uiwait(d);
    
    function pipetteResistanceAndContinue_callback(paramBtn, event)
        parametersValue = 1;
        delete(gcf)
    end
    
    function stopParameters_callback(paramBtn, event)
        parametersValue = 2;
        delete(gcf)
    end
end
