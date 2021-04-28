function parametersValue = questdlg_genericStart
    d = dialog('Position',[500, 500, 220, 60],'Name','Start Recording');

    paramBtn  = uicontrol('Parent',d,...
               'Position',[10 10 95 40],...
               'String',' REC PARAMS',...
               'Callback',@startAcquisitionParam_callback);
           
    experimBtn = uicontrol('Parent',d,...
               'Position',[115 10 40 40],...
               'String','EXIT',...
               'Callback','delete(gcf)');
    
    parametersValue = 0;       
    uiwait(d);
    
    function startAcquisitionParam_callback(paramBtn, event)
        parametersValue = 1;
        delete(gcf)
    end
    

end

