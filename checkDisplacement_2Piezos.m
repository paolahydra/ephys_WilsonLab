function allFr = checkDisplacement_2Piezos(cam, experimentHandle, varargin)

% setting schemes
pull = [-1, 1];     % (pull:   1:ANTERIOR is negative,  2:POSTERIOR is postive)
push = [1, -1];     % (push:   1:ANTERIOR is postive,  2:POSTERIOR is negative)

rampLength = 0.5; %seconds

%% set things up
HWsettings = sensor_settings_PP; % to set AO channels consistently

daqreset;
s = daq.createSession('ni');
s.Rate = 1e4;

ao = s.addAnalogOutputChannel(HWsettings.DAQdevName, ...
            [HWsettings.AOch.PIEZO_OUTPUT_LEFT, HWsettings.AOch.PIEZO_OUTPUT_RIGHT] , 'Voltage');
for ch = 1:2
    ao(ch).Range = [-10 10];
    ao(ch).TerminalConfig = 'SingleEnded';
end        


% here and generally cameras are ordered: 1 post, 2 down, 3 ant
if nargin == 3
    currentFrames = varargin{1};
    for c = 1:3
        allFr(c).frames(:,:,1) = currentFrames(c).baseline;
        allFr(c).frames(:,:,2) = currentFrames(c).piezoON;
    end
end
counter = 2;

vidnames(1).name = 'cam.post.vid';
vidnames(2).name = 'cam.down.vid';
vidnames(3).name = 'cam.ant.vid';


%% 3: baseline with piezo attached
disp('now recording frames of sample displacements. Be patient')
% capture frames
counter = counter+1;
for c = 1:3
    allFr(c).frames(:,:,counter) = getsnapshot(eval(vidnames(c).name));
end
pause(0.2)

%% 4: both push 20 um
displacement = 20;
scheme = 'push';

% new version, ramping stimulus
voltageCommand = displacement/9;
stim = linspace(0, voltageCommand, s.Rate*rampLength);
stim = stim(:);
outData = stim * eval(scheme);
queueOutputData(s, outData);
s.startForeground;
pause(0.1)

counter = counter+1;
for c = 1:3
    allFr(c).frames(:,:,counter) = getsnapshot(eval(vidnames(c).name));
end
pause(0.2)

%rampDown
queueOutputData(s, flipud(outData));
s.startForeground;
pause(0.1)

%% 5: both pull 20 um   
displacement = 20;
scheme = 'pull';

% new version, ramping stimulus
voltageCommand = displacement/9;
stim = linspace(0, voltageCommand, s.Rate*rampLength);
stim = stim(:);
outData = stim * eval(scheme);
queueOutputData(s, outData);
s.startForeground;
pause(0.1)

counter = counter+1;
for c = 1:3
    allFr(c).frames(:,:,counter) = getsnapshot(eval(vidnames(c).name));
end
pause(0.2)

%rampDown
queueOutputData(s, flipud(outData));
s.startForeground;
pause(0.1)

%% 6: baseline
counter = counter+1;
for c = 1:3
    allFr(c).frames(:,:,counter) = getsnapshot(eval(vidnames(c).name));
end
pause(0.2)

%% 7: both pull 40 um
displacement = 40;
scheme = 'pull';

% new version, ramping stimulus
voltageCommand = displacement/9;
stim = linspace(0, voltageCommand, s.Rate*rampLength*2);
stim = stim(:);
outData = stim * eval(scheme);
queueOutputData(s, outData);
s.startForeground;
pause(0.1)

counter = counter+1;
for c = 1:3
    allFr(c).frames(:,:,counter) = getsnapshot(eval(vidnames(c).name));
end
pause(0.2)

%rampDown
queueOutputData(s, flipud(outData));
s.startForeground;
pause(0.1)

%% 8: both push 40 um
displacement = 40;
scheme = 'push';

% new version, ramping stimulus
voltageCommand = displacement/9;
stim = linspace(0, voltageCommand, s.Rate*rampLength*2);
stim = stim(:);
outData = stim * eval(scheme);
queueOutputData(s, outData);
s.startForeground;
pause(0.1)

counter = counter+1;
for c = 1:3
    allFr(c).frames(:,:,counter) = getsnapshot(eval(vidnames(c).name));
end
pause(0.2)

%rampDown
queueOutputData(s, flipud(outData));
s.startForeground;
pause(0.1)

%% 9: baseline
% capture frames at steady state
counter = counter+1;
for c = 1:3
    allFr(c).frames(:,:,counter) = getsnapshot(eval(vidnames(c).name));
end
pause(0.2)


%% automatically save frames (.mat)
disp('saving timestamped frames...')
save(fullfile(experimentHandle.flyfolder, sprintf('piezosDispl_fly%d_%s.mat', experimentHandle.flyNum, datestr(now, 30))), 'allFr' )

disp('done.')

%% 
delete(s)
daqreset;

end

