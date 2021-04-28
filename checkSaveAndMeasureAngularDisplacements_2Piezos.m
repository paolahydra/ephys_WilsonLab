function [armL, armR, dispFrames] = checkSaveAndMeasureAngularDisplacements_2Piezos(cam, experimentHandle, piezosTilt_deg, allFr, varargin)
%make it more flexible
%save also the corresponding displacements and camera information
%(self-contained)
%
% input the lines
% extract the angles
% estimate arm length and output a calibration function for each of the
% piezos
% 
% 


dispFrames = [];
displacements = [0 36 0 -36 0]'; % positive is push - do not change this structure. Code is not robust.

assert(isequal(unique(displacements(1:2:end)), 0), 'code is not robust to changes in structure of variable displacements')
commands = displacements/9; % * [1, -1];  % first column is anterior, second is posterior piezo command
stimOnDur = 0.1; %seconds
rampLength = 0.9;

%% set things up

HWsettings = sensor_settings_PP; % to set AO channels consistently

daqreset;
s = daq.createSession('ni');
s.Rate = 1e4;

% ai = s.addAnalogInputChannel(HWsettings.DAQdevName, ...
%             [HWsettings.AIch.PIEZO_SENSOR_LEFT, HWsettings.AIch.PIEZO_SENSOR_RIGHT] , 'Voltage'); 
% for ch = 1:2
%     ai(ch).Range = [-10 10];
%     ai(ch).TerminalConfig = 'SingleEnded';
% end

ao = s.addAnalogOutputChannel(HWsettings.DAQdevName, ...
            [HWsettings.AOch.PIEZO_OUTPUT_LEFT, HWsettings.AOch.PIEZO_OUTPUT_RIGHT] , 'Voltage');
for ch = 1:2
    ao(ch).Range = [-10 10];
    ao(ch).TerminalConfig = 'SingleEnded';
end        


% here and generally cameras are ordered: 1 post, 2 down, 3 ant
if nargin == 5
    currentFrames = varargin{1};
    for c = 1:3
        allFr(c).frames(:,:,1) = currentFrames(c).baseline;
        allFr(c).frames(:,:,2) = currentFrames(c).piezoON;
    end
end

vidnames(1).name = 'cam.post.vid';
vidnames(2).name = 'cam.down.vid';
vidnames(3).name = 'cam.ant.vid';


%% play stimuli and record frame --  %OLD
% commands = displacements/9 * [1, -1];  % first column is anterior, second is posterior piezo command
% stimOnDur = 0.6; %seconds
% disp('stimulating and collecting frames: be patient')
% 
% % %% play stimuli and record frame -- OLD sharp stimulus
% for i = 1 : size(commands,1)    
%     disp(displacements(i))
%     voltageCommand = commands(i,:);
%     outputSingleScan(s, voltageCommand ) 
%     % capture frames at steady state
%     pause(stimOnDur)
%     for c = 1:3
%         dispFrames(c).frames(:,:,i) = getsnapshot(eval(vidnames(c).name));
%     end
%     pause(0.1)
% end

%% play stimuli and record frame - ramping
disp('stimulating and collecting frames: be patient')
% first zero
for c = 1:3
    dispFrames(c).frames(:,:,1) = getsnapshot(eval(vidnames(c).name));
end
    
for i = 2 : 2 : size(commands,1)    
    disp(displacements(i))
    voltageCommand = commands(i);
    stim = linspace(0, voltageCommand, s.Rate*rampLength);
    stim = stim(:);
    outData = stim * [1, -1];  % first column is anterior, second is posterior piezo command
    
    queueOutputData(s, outData);
    s.startForeground;
     
    % capture frames at steady state
    pause(stimOnDur)
    for c = 1:3
        dispFrames(c).frames(:,:,i) = getsnapshot(eval(vidnames(c).name));
    end
    pause(0.2)
    
    %rampDown
    queueOutputData(s, flipud(outData));
    s.startForeground;
    pause(stimOnDur)
    for c = 1:3
        dispFrames(c).frames(:,:,i+1) = getsnapshot(eval(vidnames(c).name));
    end
    pause(0.2)
end


%% play back frames and put lines on
% from the second frame on, they should take the 1st frame lines as default

numFrames = length(displacements);
posR = nan(2,2,length(displacements));
posL = nan(2,2,length(displacements));
linesR = [];
linesL = [];
colors = brewermap(length(displacements), 'Set3');


f = figure('WindowStyle', 'normal');
% f.WindowState = 'fullscreen';
f.WindowState = 'maximized';
% set(f, 'waitstatus', 'waiting')

hbox = uix.HBoxFlex( 'Parent', f ,'Spacing', 3);

    vbox1 = uix.VBox('Parent', hbox, 'Spacing', 1);
    
        handles.imSingle.ax = axes(uicontainer('Parent', vbox1),'ActivePositionProperty', 'Position');
        axis image, axis off, hold on
        
        controlsObj = uix.VButtonBox('Parent', vbox1, 'Padding',1);
        controlsObj.ButtonSize = [controlsObj.Position(3)*0.98 30];
        
            handles.sl = uicontrol('Parent',controlsObj,'Style','slider', ...
                      'value',1, 'min',1, 'max',numFrames, 'SliderStep', [1/(numFrames-1) 1/(numFrames-1)]); 
                  handles.sl.Units = 'normalized';
                  handles.sl.Position(1) = 0.02;
                  handles.sl.Position(3) = 0.96;
                  
            handles.slText = uicontrol('Parent',controlsObj, 'Style','text','String',1);

            handles.doneButton = uicontrol('Parent', controlsObj, 'Style', 'pushbutton', 'String', 'DONE');   
        
        
    vbox2 = uix.VBox('Parent', hbox, 'Spacing', 3);
        handles.imAll.ax = axes(uicontainer('Parent', vbox2),'ActivePositionProperty', 'OuterPosition');
        axis image, axis off, hold on

        hPanels = uix.HBox('Parent', vbox2, 'Spacing', 3);
            pR = uix.Panel( 'Parent', hPanels, 'Title', 'R', 'Padding', 2 );
                tR = uitable(pR,'ColumnName',{'displ';'dAng_deg';'arm_um'});
                angR_rad = nan(numFrames,1);
                armR_um = nan(numFrames,1);
                tR.Data(:, 1) = displacements;
                tR.Data(:, 2) = nan(numFrames,1);
                tR.Data(:, 3) = armR_um;
                tR.ColumnWidth = {35 67 67};
            
            pL = uix.Panel( 'Parent', hPanels, 'Title', 'L', 'Padding', 2 );
                tL = uitable(pL,'ColumnName',{'displ';'dAng_deg';'arm_um'});
                angL_rad = nan(numFrames,1);
                armL_um = nan(numFrames,1);
                tL.Data = nan(numFrames,2);
                tL.Data(:, 1) = displacements;
                tL.Data(:, 2) = nan(numFrames,1);
                tL.Data(:, 3) = armL_um;
                tL.ColumnWidth = {35 67 67};
                
       hPanels2 = uix.HBox('Parent', vbox2, 'Spacing', 3);
            handles.plotR = axes(uicontainer('Parent', hPanels2),'ActivePositionProperty', 'Position'); axis manual; hold on
            handles.plotL = axes(uicontainer('Parent', hPanels2),'ActivePositionProperty', 'Position'); axis manual; hold on
        
            
hbox.Widths = [-1.3, -1];
vbox1.Heights = [-3, -1];
vbox2.Heights = [-2.5, -1.2, -1];
linkaxes([handles.plotR, handles.plotL], 'y')
handles.plotR.YLim = [60 220];
handles.plotR.XLim = [min(displacements), max(displacements)];
handles.plotL.XLim = [min(displacements), max(displacements)];


%% set position 1 as default 
markedFrames = false(1, numFrames);

posR(:,:,1) = [ 530  500;  225  640];  %change default to something more plausible
posL(:,:,1) = [1060  500; 1370  640];

axes(handles.imSingle.ax)
imshow(dispFrames(2).frames(:,:,1), [])

lR = imline(gca, posR(:,:,1)); 
lL = imline(gca, posL(:,:,1)); %base x-y; tip x-y
setColor(lR, 'y');
setColor(lL, 'y');
previousFrame = 1;

axes(handles.imAll.ax)
imshow(dispFrames(2).frames(:,:,1), [])


lh  = addlistener(handles.sl,'Action',@(src,evnt)surfzlim(src,evnt,handles));
lh1 = addlistener(handles.doneButton, 'Action', @(src,evnt)doneOut(src,evnt,handles));
        
uiwait(f);

%%
function surfzlim(source, event, handles)
    val = floor(source.Value);
    handles.slText.String = num2str(val);

    
    % first, update/store previous frame position info 
    posR(:,:,previousFrame) = lR.getPosition;
    posL(:,:,previousFrame) = lL.getPosition;
    lR.delete;
    lL.delete;
    markedFrames(previousFrame) = true;
    
    axes(handles.imAll.ax)
    linesR(previousFrame).l = line(posR(:,1,previousFrame), posR(:,2,previousFrame), 'Color', colors(previousFrame,:));
    linesL(previousFrame).l = line(posL(:,1,previousFrame), posL(:,2,previousFrame), 'Color', colors(previousFrame,:));
    linesR(previousFrame).l.DisplayName = sprintf('%2d',previousFrame);
    linesL(previousFrame).l.DisplayName = sprintf('%2d',previousFrame);
    
    grouping = cat(2, find(markedFrames), find(markedFrames));
    clickableLegend(    cat(1,linesR.l, linesL.l),    'groups', grouping, 'Location', 'northeastoutside');
    
    
    %% 
    baselines = displacements == 0;
    ibaselines = find(baselines);
    
    ang_rad_zero = [];
    for b = 1 : sum(squeeze(~isnan(posR(1,1,baselines))))
        b_i = ibaselines(b);
%         ang_rad_zero(b) = atan((posR(2,2,b_i)-posR(1,2,b_i)) / (posR(2,1,b_i)-posR(1,1,b_i)));
        ang_rad_zero(b) = atan((posR(2,2,b_i)-posR(1,2,b_i)) / (posR(2,1,b_i)-posR(1,1,b_i))) / cosd(piezosTilt_deg.R);
    end
    angR_rad_zero = circ_mean(ang_rad_zero(:));

    
    ang_rad_zero = [];
    for b = 1 : sum(squeeze(~isnan(posL(1,1,baselines))))
        b_i = ibaselines(b);
%         ang_rad_zero(b) = atan((posL(2,2,b_i)-posL(1,2,b_i)) / (posL(2,1,b_i)-posL(1,1,b_i)));
        ang_rad_zero(b) = atan((posL(2,2,b_i)-posL(1,2,b_i)) / (posL(2,1,b_i)-posL(1,1,b_i))) / cosd(piezosTilt_deg.L);
    end
    angL_rad_zero = - circ_mean(ang_rad_zero(:)); % note sign for L
    
    
    % for all the frames, resubtract averaged baseline
    % for frames w
    for b = 1 : length(displacements)
%         angR_rad(b,1) = atan((posR(2,2,b)-posR(1,2,b)) / (posR(2,1,b)-posR(1,1,b))) - angR_rad_zero;    
%         angL_rad(b,1) = - atan((posL(2,2,b)-posL(1,2,b)) / (posL(2,1,b)-posL(1,1,b))) - angL_rad_zero; % note sign for L
        angR_rad(b,1) = atan((posR(2,2,b)-posR(1,2,b)) / (posR(2,1,b)-posR(1,1,b))) / cosd(piezosTilt_deg.R) - angR_rad_zero;    
        angL_rad(b,1) = - atan((posL(2,2,b)-posL(1,2,b)) / (posL(2,1,b)-posL(1,1,b))) / cosd(piezosTilt_deg.L) - angL_rad_zero; % note sign for L
        
        % if not a zero displacement frame
        if displacements(b) ~= 0
            armR_um(b) = displacements(b)/tan(angR_rad(b)); %fix others too
            armL_um(b) = displacements(b)/tan(angL_rad(b));
        end
    end
    
    % always update the entire table, in case user updates frame 1 later on
    angR = rad2deg(angR_rad);
    angL = rad2deg(angL_rad);
    tR.Data(:, 2) = round(angR*100)/100;
    tR.Data(:, 3) = armR_um;
    
    tL.Data(:, 2) = round(angL*100)/100;
    tL.Data(:, 3) = armL_um;
    
    
    %% then update axes and process new frame
    axes(handles.imSingle.ax)
    
    handles.imSingle.ax.NextPlot = 'replaceChildren';
    imshow(dispFrames(2).frames(:,:,val), [])
    handles.imSingle.ax.NextPlot = 'add';
    
    % display defaults
    if isnan(posR(1,1,val))
        lR = imline(gca, posR(:,:,1));
        lL = imline(gca, posL(:,:,1));
        setColor(lR, 'y');
        setColor(lL, 'y');
        
    elseif ~isnan(posR(1,1,val)) % pre-existing line 
        lR = imline(gca, posR(:,:,val));
        lL = imline(gca, posL(:,:,val));
        setColor(lR, 'g');
        setColor(lL, 'g');

    end
    
    
    axes(handles.plotR)
    handles.plotR.NextPlot = 'replaceChildren';
    plot(displacements, armR_um, 'xb', 'MarkerSize',12)
    axes(handles.plotL)
    handles.plotL.NextPlot = 'replaceChildren';
    plot(displacements, armL_um, 'xr', 'MarkerSize',12)
    
    % before you leave this function
    previousFrame = val;
end

%%
function doneOut(source, event, handles)
    
    % save last line
    posR(:,:,previousFrame) = lR.getPosition;
    posL(:,:,previousFrame) = lL.getPosition;
    lR.delete;
    lL.delete;
    markedFrames(previousFrame) = true;
    
    axes(handles.imAll.ax)
    linesR(previousFrame).l = line(posR(:,1,previousFrame), posR(:,2,previousFrame), 'Color', colors(previousFrame,:));
    linesL(previousFrame).l = line(posL(:,1,previousFrame), posL(:,2,previousFrame), 'Color', colors(previousFrame,:));
    linesR(previousFrame).l.DisplayName = sprintf('%2d',previousFrame);
    linesL(previousFrame).l.DisplayName = sprintf('%2d',previousFrame);
   
    
    %% recalculate angles zero from all baselines available, and all angles in Table
    baselines = displacements == 0;
    ibaselines = find(baselines);
    
    ang_rad_zero = [];
    for b = 1 : sum(squeeze(~isnan(posR(1,1,baselines))))
        b_i = ibaselines(b);
%         ang_rad_zero(b) = atan((posR(2,2,b_i)-posR(1,2,b_i)) / (posR(2,1,b_i)-posR(1,1,b_i)));
       ang_rad_zero(b) = atan((posR(2,2,b_i)-posR(1,2,b_i)) / (posR(2,1,b_i)-posR(1,1,b_i))) / cosd(piezosTilt_deg.R);
    end
    angR_rad_zero = circ_mean(ang_rad_zero(:));

    
    ang_rad_zero = [];
    for b = 1 : sum(squeeze(~isnan(posL(1,1,baselines))))
        b_i = ibaselines(b);
%         ang_rad_zero(b) = atan((posL(2,2,b_i)-posL(1,2,b_i)) / (posL(2,1,b_i)-posL(1,1,b_i)));
        ang_rad_zero(b) = atan((posL(2,2,b_i)-posL(1,2,b_i)) / (posL(2,1,b_i)-posL(1,1,b_i))) / cosd(piezosTilt_deg.L);
    end
    angL_rad_zero = - circ_mean(ang_rad_zero(:)); % note sign for L
    
    
    % for all the frames, resubtract averaged baseline
    % for frames w
    for b = 1 : length(displacements)
%         angR_rad(b,1) = atan((posR(2,2,b)-posR(1,2,b)) / (posR(2,1,b)-posR(1,1,b)))  - angR_rad_zero;    
%         angL_rad(b,1) = - atan((posL(2,2,b)-posL(1,2,b)) / (posL(2,1,b)-posL(1,1,b))) - angL_rad_zero; % note sign for L
        angR_rad(b,1) = atan((posR(2,2,b)-posR(1,2,b)) / (posR(2,1,b)-posR(1,1,b))) / cosd(piezosTilt_deg.R) - angR_rad_zero;    
        angL_rad(b,1) = - atan((posL(2,2,b)-posL(1,2,b)) / (posL(2,1,b)-posL(1,1,b))) / cosd(piezosTilt_deg.L) - angL_rad_zero; % note sign for L

        % if not a zero displacement frame
        if displacements(b) ~= 0
            armR_um(b) = displacements(b)/tan(angR_rad(b)); %fix others too
            armL_um(b) = displacements(b)/tan(angL_rad(b));
        end
    end
    
    % always update the entire table, in case user updates frame 1 later on
    angR = rad2deg(angR_rad);
    angL = rad2deg(angL_rad);
    tR.Data(:, 2) = round(angR*100)/100;
    tR.Data(:, 3) = armR_um;
    
    tL.Data(:, 2) = round(angL*100)/100;
    tL.Data(:, 3) = armL_um;
    
    
    armL = nanmedian(armL_um); % median
    armR = nanmedian(armR_um); % median
    
    % save all 
    disp('saving timestamped frames and data...')
    save(fullfile(experimentHandle.flyfolder, sprintf('piezosDispl_fly%d_%s.mat', experimentHandle.flyNum, datestr(now, 30))), 'posR', 'posL', 'armR_um', 'armL_um', 'angR_rad', 'angL_rad',   ...
                                                                                                                              'dispFrames', 'displacements', 'commands', 'stimOnDur', 'allFr' )
    if isequal(get(f, 'waitstatus'), 'waiting')
        %the GUI is still in UIWAIT, use UIRESUME
        uiresume(f);
    end
    delete(s)
    daqreset; 
    delete(f)
end

end

