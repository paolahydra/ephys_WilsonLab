function [stopValue, d] = questdlg_PPStop
    d = dialog('Position',[2951, 1080, 100, 50],'Name','Start Parameters Acquisition');
    
    stopBtn = uicontrol('Parent',d,...
               'Position',[25 12 70 25],...
               'String','STOP',...
               'Callback','delete(gcf)');
    
    stopValue = 1;       
    uiwait(d);
    
end

