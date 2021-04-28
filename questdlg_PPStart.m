function parametersValue = questdlg_PPStart
    d = dialog('Position',[2951, 1080, 270, 60],'Name','Start Recording'); % [2951, 1520, 270, 60]

    paramBtn  = uicontrol('Parent',d,...
               'Position',[10 10 95 40],...
               'String',' REC PARAMS',...
               'Callback',@startAcquisitionParam_callback);
           
    experimBtn = uicontrol('Parent',d,...
               'Position',[115 10 95 40],...
               'String','REC EXPERIM',...
               'Callback',@startEperiment_callback);
           
    experimBtn = uicontrol('Parent',d,...
               'Position',[220 10 40 40],...
               'String','EXIT',...
               'Callback','delete(gcf)');
    
    parametersValue = 0;       
    uiwait(d);
    
    function startAcquisitionParam_callback(paramBtn, event)
        parametersValue = 1;
        delete(gcf)
    end
    
    function startEperiment_callback(paramBtn, event)
        parametersValue = 2;
        delete(gcf)
    end
end

