%% [BEFORE-FLY] experiment setup
experimentHandle = getNewFlyID('D:\Dropbox (HMS)\p2');  % lets you choose between last run, new run or new fly
load('D:\Dropbox (HMS)\p2\antennalMovements\flies4to9\modelsAngularData.mat', 'angularData');
armR            = 210; %default
armL            = 210; %default
pipetteNum      = 0; %start from here, will be updated every time parameters are recorded
DAQ_fs          = 4e4;
metadata.key.ID = datestr(now, 30); %default, changed when new metadata are saved

% experimentHandle = getNewFlyID(rootfolder, lastFly)   % optional syntax to override last fly specification
%% [BEFORE-FLY] set cameras up
% cameras are ordered: 1 post, 2 down, 3 ant
imaqreset
set(0, 'DefaultFigureWindowStyle', 'docked')
cam = setupcameras_experiment();

%% [BEFORE-FLY P2] set piezo DC-offset
DC_offset = setDCoffset_2Piezos;     % resets and releases any DAC that might be connected 



%click!




%% patch before attaching piezos
run_continuousAcquisition;


%% optional - NEW CELL
experimentHandle = getNewCellID(experimentHandle);  % lets you choose between last cell, new cell or new fly






%% position piezos 
assert(logical(exist('DC_offset', 'var')), 'set piezos DC offset first')
disp('Use GUI:')
currentFrames2 = positionPiezos(cam); %check for GUI - o saving here



%% visualize piezo displacement / antennal rotation
% all frames are automatically saved with timestamps in flyfolder
if exist('currentFrames2', 'var')
    allFr = checkDisplacement_2Piezos(cam, experimentHandle, currentFrames2); % 12, 24 um bidirectional
else
    allFr = checkDisplacement_2Piezos(cam, experimentHandle); % 12, 24 um bidirectional
end
% ok, now visualize them: for now this will do:
warning('off')
figure; imshowpair(allFr(2).frames(:,:,7), allFr(2).frames(:,:,8)) %biggest push/pull displacements - change counts...
figure; imshowpair(allFr(1).frames(:,:,7), allFr(1).frames(:,:,8)) %biggest push/pull displacements - change counts...
figure; imshowpair(allFr(3).frames(:,:,7), allFr(3).frames(:,:,8)) %biggest push/pull displacements - change counts...
warning('on')



%% calibrate stimuli bilaterally
piezosTilt_deg.R = 25; %this should really be the antenna's tilt, not piezo's - but they are overall fairly matched
piezosTilt_deg.L = 17; %currently used - this was a good call

% now using bigger displacements too (36um vs 28um)
[armLor, armRor] = checkSaveAndMeasureAngularDisplacements_2Piezos(cam, experimentHandle, piezosTilt_deg, allFr); %DOUBLECHECK CORRECTIONS HERE
fprintf('L arm: %4.0f um \t - L max displ - originally: %2.1f um\n', armLor, armLor * tand(12.4) ) %12.4 is the biggest angular displacement in the stim set, approximately...
fprintf('R arm: %4.0f um \t - R max displ - originally: %2.1f um\n', armRor, armRor * tand(12.4) )

% % in use up to 2/4
m = max([armRor, armLor]);
% fact2 = 44.97 / (m*tand(12.4));
% fact = min([fact2, 1.2]);
% 
% % fact = 1.02;
% armR = fact * armRor;
% armL = fact * armLor; % I would need to make the stimuli even bigger, 
% %    because the measured ant movements tuning curves were an underestimate 
% %    of actual antennal rotation (never corrected by the tilt angle)
% % however, I am hitting the piezo range already...........

% target shorter arms if possible. keep em under 193 um
if m > 193
    error('reduce the arm')
end



disp('-----------------------')
fprintf('L max displ - now: %2.1f um\n', armL * tand(12.4) ) % 12.3772 is the biggest angular displacem used
fprintf('R max displ - now: %2.1f um\n', armR * tand(12.4) )

%% [POST-PATCHING CELL X] prepare piezo stims - all automatically saved
% pipetteNum = 0; %start from here, will be updated every time parameters are recorded
% DAQ_fs              = 4e4;
% % assert(armL * tand(max(abs(ang_displacements))) <= 44 && armR * tand(max(abs(ang_displacements))) <= 44, ...
% %     'stimuli out of piezo limits. Move the probe closer to the midline')

N_withinPseudoBlock = 1;    % 2;
N_nPseudoBlocks     = 5;    % 8;    % with the ability to stop earlier if need so
singleStimDur       = 1.3;  % 1.25; % seconds

factor_jointZero    = 4;    % integer. Set to 1 to have the same number of baselines as all the other stimuli. 2 for two as many, and so forth.

% angDownsampling = 3; %3 keeps every 15 degrees
% iWind = angularData.fly5R.i30;
% metadata = bilateralPiezoCommands_windPlayback_ang(iWind, angDownsampling, N_withinPseudoBlock,N_nPseudoBlocks,factor_jointZero,singleStimDur, DAQ_fs, experimentHandle, armL, armR); 
% iWind = angularData.fly9;

metadata = bilateralPiezoCommands_windPlaybackFULL_ang(N_withinPseudoBlock, factor_jointZero, N_nPseudoBlocks,singleStimDur, DAQ_fs, experimentHandle, armL, armR);


% run experiment (click)
run_continuousAcquisition;


% click!
 

%% [DONE FOR THE DAY] close all cameras and exit (save!!!)
stopClosePreviews(cam)


%% [DONE FOR THE DAY] exit (and save scripts when asked)
exit







%% [post] read in parameters data file
HWsettings = sensor_settings_PP; % to set AO channels consistently
[filename, pathname] = uigetfile('parameters_*.bin', 'Select file of parameters to load');
fid = fopen(fullfile(pathname, filename), 'r');
data = fread(fid, 'double' );
fclose(fid);

data = reshape(data, 7, []);
timestamps = (data(1,:))';
data = data(2:end,:);
data = data';  %for consistency with TO's code

figure; plot(timestamps(1:4001), data( 1:4001, HWsettings.DATAch.VOLTAGE_NON_SCALED_A_DAQ_AI ))
hold on

[ current, voltage, scaled ] = get_scaled_voltage_and_current_PP( data );
figure; hold on
plot(timestamps(1:4001), scaled(1:4001)); ylabel('scaled, [pA]')
plot(timestamps(1:4001), current(1:4001));
% 
% 
% figure; hold on
% plot(timestamps, scaled); ylabel('scaled, [pA]')
% plot(timestamps, current);
