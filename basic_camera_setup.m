%% prelims
%  install  point grey support package from the Image Acquisition Toolbox
% https://www.mathworks.com/hardware-support/point-grey-camera.html


%  matlabcentral fileexchange used:
% https://www.mathworks.com/matlabcentral/fileexchange/31437-windowapi


%  ImageJ plugin to open avi files:
% https://imagej.nih.gov/ij/plugins/file-handler.html


%  ImageAcqToolbox tutorial:
% https://www.mathworks.com/help/imaq/basic-image-acquisition-procedure.html
% 
%  see Programmatic Workflow:
% https://www.mathworks.com/help/imaq/acquisition-using-any-hardware.html


%  PP - 3/15/2018

 
%% get cameras
imaqhwinfo 

% adaptorName = 'gige';       % does not work
adaptorName = 'pointgrey';    % 'pointgrey' sometimes only sees two cameras (but now it sees all 4)
hwInfo = imaqhwinfo(adaptorName); 
cameraNames = cat(2, hwInfo.DeviceInfo.DeviceName);
cameraNums = length(hwInfo.DeviceInfo);
cameraNamesReadable = cat(1, hwInfo.DeviceInfo.DeviceName);
disp(cameraNamesReadable)

%% create a video input object and set properties
pause(1)
ID_camera = ceil(strfind(cameraNames, '12A2M') / (length(cameraNames)/cameraNums)); %name will still be ambiguous for same model cameras.
pause(1)
% hwInfo.DeviceInfo(ID_camera)

vid = videoinput(adaptorName, ID_camera, hwInfo.DeviceInfo(ID_camera).DefaultFormat);
% get(vid)        % The get function lists all the properties of the object with their current values.

src = getselectedsource(vid);
% src.SerialNumber % use it to disambiguate the two 13X cameras


%% Preview the Video Stream (at any time after you create the vid object)   
handles.Image = preview(vid);
handles.Figure = ancestor(handles.Image, 'figure');
handles.Figure.WindowStyle = 'normal';

handles.Figure.Units = 'pixels';
WindowAPI(handles.Figure,'Position', 'full', 1);  % to make it a *real* fullscreen, in the monitor specified by last argument
handles.Figure.Units = 'pixels'; 
WindowAPI(handles.Figure,'Position', 'full', 1);  % at least one of these two lines does not work the first time...

% handles.Axes = ancestor(handles.Image, 'axes');
% % this so far works but the magnification is fixed to 1 and cannot be changed.
% % If the framesize is bigger than your monitor, you can at least specify your desired vid.ROIPosition
% % 
% % Alternatively, this:
% handles.Axes.Units = 'normalized';
% handles.Axes.Position = [0 0 1 1]; 
% % would work to stretch the image to fit the screen, but would have to be
% % called every time the figure is updated (online) in the preview
% % function, and I have not tried to implement it.
% settings

% not sure why, but I have to do it after starting the preview
src.FrameRate = 19;
src.Gamma = 1;



%% Stop the preview
stoppreview(vid)

% % optional:
closepreview(vid)   
% closepreview          % Calling closepreview without any arguments closes all open Video Preview windows.



%% acquire a single frame to workspace
% single frame into matlab workspace
frame = getsnapshot(vid);

% other acquisition examples:

% %% acquire multiple frames, manually triggering
% % triggerinfo(vid)      % get possible trigger configurations
% triggerconfig(vid, 'Manual')
% vid.FramesPerTrigger = 12;
% 
% start(vid), wait(vid);
% trigger(vid)
% flushdata(vid,'triggers'); %delete data buffered before the trigger.
% 
% %% acquire multiple frames, starting immediately
% % configs = triggerinfo(vid)      % get possible trigger configurations
% triggerconfig(vid, 'immediate')
% vid.FramesPerTrigger = 14;
% 
% start(vid), wait(vid);
% 
% 
% %% check status (optional) and retrieve data from memory
% 
% isrunning(vid)
% islogging(vid)
% 
% data = getdata(vid); % alternatively, flushdata(vid)
% size(data)
% 
% %% where do you want to save?
% vid.LoggingMode = 'memory'; % default, faster but limited.
% vid.LoggingMode = 'disk';   % https://www.mathworks.com/help/imaq/logging-image-data-to-disk.html -- READ
% vid.LoggingMode = 'disk&memory';




%% save multiple frames to disk, starting immediately
% https://www.mathworks.com/help/imaq/logging-image-data-to-disk.html

% configs = triggerinfo(vid)      % get possible trigger configurations
triggerconfig(vid, 'immediate')
vid.LoggingMode = 'disk';
logfile = VideoWriter('logfile.avi'); %'Motion JPEG AVI' default compression

vid.DiskLogger = logfile; % do not modify the VideoWriter object past this point

% When the FramesPerTrigger property is set to Inf, the object ignores the value of the TriggerRepeat property
vid.FramesPerTrigger = inf; % acquire until stop, otherwise specify a number
vid.FrameGrabInterval = 2; %if you want to save at lower fr than previewing

start(vid) %this overwrites previously saved data at the same logfile if a new one is not supplied

pause(5) % or do other stuff (e.g. stimulating antennae)

stop(vid) % when done

while (vid.FramesAcquired ~= vid.DiskLoggerFrameCount) 
    pause(.1)
end



%% camera 3
ID_camera3 = ceil(strfind(cameraNames, '20E4M') / (length(cameraNames)/cameraNums)); %name will still be ambiguous for same model cameras.
% hwInfo.DeviceInfo(ID_camera3)

vid3 = videoinput(adaptorName, ID_camera3, hwInfo.DeviceInfo(ID_camera3).DefaultFormat);
% get(vid3)        % The get function lists all the properties of the object with their current values.

src = getselectedsource(vid3);
% src.SerialNumber % use it to disambiguate the two 13X cameras


% Preview the Video Stream (at any time after you create the vid object)   
handles.Image = preview(vid3);
handles.Figure = ancestor(handles.Image, 'figure');
handles.Figure.WindowStyle = 'normal';

handles.Figure.Units = 'pixels';
WindowAPI(handles.Figure,'Position', 'work', 3);  % to make it a *real* fullscreen, in the monitor specified by last argument
handles.Figure.Units = 'pixels'; 
WindowAPI(handles.Figure,'Position', 'work', 3);  % at least one of these two lines does not work the first time...


% stoppreview(vid3)
% closepreview(vid3)  



%% camera Right antenna (posterior piezo)
ID_camera2 = ceil(strfind(cameraNames, '13E4M') / (length(cameraNames)/cameraNums)); %name will still be ambiguous for same model cameras.
vid2 = videoinput(adaptorName, ID_camera2(1), hwInfo.DeviceInfo(ID_camera2(1)).DefaultFormat);

src = getselectedsource(vid2);
SerialNumber = src.SerialNumber;
delete src
delete(vid2)
clear vid2 src

if strcmp(SerialNumber, '14466065') % use it to disambiguate the two 13X cameras
    vidPosterior = videoinput(adaptorName, ID_camera2(1), hwInfo.DeviceInfo(ID_camera2(1)).DefaultFormat);
    vidAnterior  = videoinput(adaptorName, ID_camera2(2), hwInfo.DeviceInfo(ID_camera2(2)).DefaultFormat);
else
    vidAnterior = videoinput(adaptorName, ID_camera2(1), hwInfo.DeviceInfo(ID_camera2(1)).DefaultFormat);
    vidPosterior  = videoinput(adaptorName, ID_camera2(2), hwInfo.DeviceInfo(ID_camera2(2)).DefaultFormat);
end
    
% Preview the Video Stream (at any time after you create the vid object)   
handlesP.Image = preview(vidPosterior);
handlesP.Figure = ancestor(handlesP.Image, 'figure');
% handlesP.Figure.WindowStyle = 'normal';
% handlesP.Figure.Units = 'pixels';
% WindowAPI(handlesP.Figure,'Position', 'work', 4);  % to make it a *real* fullscreen, in the monitor specified by last argument
% handlesP.Figure.Units = 'pixels'; 
% WindowAPI(handlesP.Figure,'Position', 'work', 4);  % at least one of these two lines does not work the first time...
srcP = getselectedsource(vidPosterior);
srcP.FrameRate = 40;
srcP.Gamma = 1;
srcP.Exposure = 0.3;
pause(0.2)
srcP.FrameRate = 40;
srcP.Gamma = 1;
srcP.Exposure = 0.3;


handlesA.Image = preview(vidAnterior);
handlesA.Figure = ancestor(handlesA.Image, 'figure');
% handlesA.Figure.WindowStyle = 'normal';
% handlesA.Figure.Units = 'pixels';
% WindowAPI(handlesA.Figure,'Position', 'work', 4);  % to make it a *real* fullscreen, in the monitor specified by last argument
% handlesA.Figure.Units = 'pixels'; 
% WindowAPI(handlesA.Figure,'Position', 'work', 4);  % at least one of these two lines does not work the first time...
srcA = getselectedsource(vidAnterior);
srcA.FrameRate = 40;
srcA.Gamma = 1;
srcA.Exposure = 0.3;
pause(0.2)
srcA.FrameRate = 40;
srcA.Gamma = 1;
srcA.Exposure = 0.3;


%% grab a frame each
frame_Up = getsnapshot(vid);
frame_Down = getsnapshot(vid3);
frame_Anterior = getsnapshot(vidAnterior);
frame_Posterior = getsnapshot(vidPosterior);


%% concatenate a frame each
frame_Up = cat(3, frame_Up, getsnapshot(vid));
frame_Down = cat(3, frame_Down, getsnapshot(vid3));
frame_Anterior = cat(3, frame_Anterior, getsnapshot(vidAnterior));
frame_Posterior = cat(3, frame_Posterior, getsnapshot(vidPosterior));

%% save all frames in current folder
save(sprintf('frames_%s.mat',datestr(now, 30)), 'frame_Up', 'frame_Down', 'frame_Anterior', 'frame_Posterior' )

%% clean up (always)
stoppreview(vid)
closepreview(vid)  

stoppreview(vid3)
closepreview(vid3)  

stoppreview(vidPosterior)
closepreview(vidPosterior)  

stoppreview(vidAnterior)
closepreview(vidAnterior)  


delete(vid)
delete(vidAnterior)
delete(vidPosterior)
delete(vid3)

delete(src) 
delete(srcA)
delete(srcP)

clear vid vidAnterior vidPosterior vid3 src srcA srcP