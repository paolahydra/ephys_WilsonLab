function [stimset] = bilateralPiezoCommands_windPlaybackFULL_ang(N_withinPseudoBlock, factor_jointZero, N_nPseudoBlocks, singleStimDur, sample_fs, experimentHandle, armL, armR )
% 190205 - corrected indexing issue that assigned rest to stim numbered 1 2
% 3 and 4. OK_PP
% 190207 - introduced a cosine-correction for big angle/displacement ratio:
% divide command of displacemene by the cosine of the angle to produce
% bigger displacements for intented bigger angles.
%
%
% negative displacements is always pull!
%
% inputs:
%
% iWind has two fields (L, R) and corresponds to one specific wind intensity of
% a fly model (select it and feed it to the function). iWind is sampled
% every 5 degrees (length of each L/R vector: 61)

% angDownsampling, e.g.: 3, to sample every 15 degrees only

% N_withinPseudoBlock = 4;
% factor_jointZero = 2;     % integer. Set to 1 to have the same number of baselines as all the other stimuli. 2 for two as many, and so forth.
% N_nPseudoBlocks = 2;      % with the ability to stop earlier if need so
% singleStimDur = 1.25;     % seconds
% sample_fs = 4e4; 
%

load('D:\Dropbox (HMS)\p2\antennalMovements\flies4to9\modelsAngularData.mat', 'samples', 'angularData');
plotfigures = 0;

% I need the model, and the x and y intersects of the 30 60 and 120 intensities 
% (generally: fly 9 model)

% fly 9:
iWind = angularData.fly9;
cf30 = samples(6,1).cf;
cf60 = samples(6,2).cf;
cf120 = samples(6,4).cf;

% % fly 8
% iWind = angularData.fly8L;
% cf30 = samples(5,1).cf_Lonly;
% cf60 = samples(5,2).cf_Lonly;
% cf120 = samples(5,4).cf_Lonly;
% 
% % fly 5
% iWind = angularData.fly5R;
% cf30 = samples(2,1).cf_Ronly;
% cf60 = samples(2,2).cf_Ronly;
% cf120 = samples(2,4).cf_Ronly;





% %% build type 2018/12/05
% 
% % build stim set up
% angL = [];
% angR = [];
% directions = [];
% intensities = [];
% 
% windDirection = -150:5:150;
% assert(length(windDirection) == length(iWind.i120.L))
% downSamplingVector = 1: angDownsampling: length(windDirection);
% 
% 
% % 120 stims
% windIntensity = 120;
% angL_build = iWind.i120.L(downSamplingVector);
% angR_build = iWind.i120.R(downSamplingVector);
% directions_build = windDirection(downSamplingVector)';
% intensities_build = windIntensity*ones(length(downSamplingVector),1);
% 
% angL = cat(1, angL, angL_build(:));
% angR = cat(1, angR, angR_build(:));
% directions = cat(1, directions, directions_build);
% intensities = cat(1, intensities, intensities_build);
% 
% % 60 stims
% windIntensity = 60;
% angL_build = iWind.i60.L(downSamplingVector);
% angR_build = iWind.i60.R(downSamplingVector);
% directions_build = windDirection(downSamplingVector)';
% intensities_build = windIntensity*ones(length(downSamplingVector),1);
% 
% angL = cat(1, angL, angL_build(:));
% angR = cat(1, angR, angR_build(:));
% directions = cat(1, directions, directions_build);
% intensities = cat(1, intensities, intensities_build);
% 
% % 30 stims
% downSamplingVector = find(ismember(windDirection, [-90:15:-60, -30:15:30, 60:15:90]));
% windIntensity = 30;
% angL_build = iWind.i30.L(downSamplingVector);
% angR_build = iWind.i30.R(downSamplingVector);
% directions_build = windDirection(downSamplingVector)';
% intensities_build = windIntensity*ones(length(downSamplingVector),1);
% 
% angL = cat(1, angL, angL_build(:));
% angR = cat(1, angR, angR_build(:));
% directions = cat(1, directions, directions_build);
% intensities = cat(1, intensities, intensities_build);
% 
% 
% %include monoantennal stimulations, doubling the zero stimulus
% xR = -150:0.001:150;
% aL30 = flipud(cf30(xR));
% aL60 = flipud(cf60(xR));
% aL120 = flipud(cf120(xR));
%     
%     
% %R antenna at zeros
% %wind directions (x points)
% wds = [fzero(cf120, -70), fzero(cf60, -70), fzero(cf30, -70), ... 
%        fzero(cf30, 70), fzero(cf60, 70), fzero(cf120, 70)]; %add [0,0]
% % find corresponding l antennal movements (at the same x points)
% x_ind = [];
% for i = 1:length(wds)
%     x_ind(i) = find((xR - wds(i)) > 0, 1);
% end
% angL_build = [aL120(x_ind(1)), aL60(x_ind(2)), aL30(x_ind(3)),...
%               aL30(x_ind(4)), aL60(x_ind(5)), aL120(x_ind(6)) ]';
% angR_build = zeros(size(angL_build));
% directions_build = wds(:);
% intensities_build = [120 60 30 30 60 120]';
% 
% angL = cat(1, angL, angL_build(:));
% angR = cat(1, angR, angR_build(:));
% directions = cat(1, directions, directions_build);
% intensities = cat(1, intensities, intensities_build);
% 
% 
% %L antenna at zeros
% %wind directions (x points)
% angL_build = zeros(size(angL_build));
% angR_build = [aL120(x_ind(1)), aL60(x_ind(2)), aL30(x_ind(3)),...
%               aL30(x_ind(4)), aL60(x_ind(5)), aL120(x_ind(6)) ]';
% 
% angL = cat(1, angL, angL_build(:));
% angR = cat(1, angR, angR_build(:));
% directions = cat(1, directions, directions_build);
% intensities = cat(1, intensities, intensities_build);
% 
% % add (0,0) now, factored, at the beginning
% angL = cat(1, zeros(factor_jointZero,1), angL);
% angR = cat(1, zeros(factor_jointZero,1), angR);
% directions = cat(1, 999*ones(factor_jointZero,1), directions);
% intensities = cat(1, zeros(factor_jointZero,1), intensities);
% 
% 
% 
% 
% displacementsL = armL * tand(angL);
% displacementsR = armR * tand(angR);



%% try different build type 2018/12/11

% build stim set up
angL = [];
angR = [];
directions = [];
intensities = [];


% new 60-isoI sample points
windIntensity = 60;
dirs = sort([0,15,22.5,30,39.5, 45, 60, 65, 70:10:100, 135, -15,-22.5,-30,-39.5, -45, -60:-10:-90, -105, -135]); 
angR_build = cf60(dirs); %yes. flip sign of directions and sort accordingly for L
angL_build = (cf60(-dirs));
directions_build = dirs';
intensities_build = windIntensity*ones(length(dirs),1);

angL = cat(1, angL, angL_build(:));
angR = cat(1, angR, angR_build(:));
directions = cat(1, directions, directions_build);
intensities = cat(1, intensities, intensities_build);


% 120 stims
windIntensity = 120;
dirs = sort([0,22.5,30,39.5, 45, 48.25,  60, 65, 70:5:100, 135, -22.5,-30,-39.5, -45, -48.25, -55:-5:-90, -115, -135, -150]); %no 15, added 48
%points are generally shifted with respect to true linear projections from
%60I, but in the right direction to make stronger claims, if true effect, or to
%leave it to noise if uncertain.
angR_build = cf120(dirs); %yes. flip sign of directions and sort accordingly for L
angL_build = (cf120(-dirs));
directions_build = dirs';
intensities_build = windIntensity*ones(length(dirs),1);

angL = cat(1, angL, angL_build(:));
angR = cat(1, angR, angR_build(:));
directions = cat(1, directions, directions_build);
intensities = cat(1, intensities, intensities_build);


% 30 stims
windIntensity = 30;
dirs = sort([15:15:30, 39.5, 60:15:90, -90, -74.6, -39.5, -30, -15]);
angR_build = cf30(dirs); %yes. flip sign of directions and sort accordingly for L
angL_build = (cf30(-dirs));
directions_build = dirs';
intensities_build = windIntensity*ones(length(dirs),1);

angL = cat(1, angL, angL_build(:));
angR = cat(1, angR, angR_build(:));
directions = cat(1, directions, directions_build);
intensities = cat(1, intensities, intensities_build);

if plotfigures
    figure; scatter(angL, angR, 'filled'); axis image; hold on
    length(directions)
    hold on
%     plot([-10 10], [-10 10])
%     plot([-10 10], [10 -10])
end

%include monoantennal stimulations, doubling the zero stimulus --
%DIFFERENT - directly compare to conjunctions
% singZeros = sort([-cf120(105), -cf60(115.55), cf30(-48.5), cf60(115.55)]); %first negative is taken from cf120(105), included already
% singZeros = sort([-cf120(0), cf60(115.55), cf30(121.4), cf30(-48.5), -cf60(115.55), cf120(0)]);
singZeros = sort([-cf120(0), -cf60(0), cf30(121.4), cf30(-48.5),  cf60(0),      cf120(0)]);


% move Left antenna only
angL_build = singZeros';
angR_build = zeros(size(singZeros));
directions_build = 999*ones(length(singZeros), 1);
intensities_build = [-3:-1, 1:3]; % fixed 19/01/07

angL = cat(1, angL, angL_build(:));
angR = cat(1, angR, angR_build(:));
directions = cat(1, directions, directions_build);
intensities = cat(1, intensities, intensities_build(:));

if plotfigures
    scatter(angL_build, angR_build, 'xr');
        length(directions)
    hold on
%     plot([-10 10], [-10 10])
%     plot([-10 10], [10 -10])
end

% move Right antenna only
angL_build = zeros(size(singZeros));
angR_build = singZeros';
directions_build = 999*ones(length(singZeros), 1);
intensities_build = [-3:-1, 1:3]; % fixed 19/01/07

angL = cat(1, angL, angL_build(:));
angR = cat(1, angR, angR_build(:));
directions = cat(1, directions, directions_build);
intensities = cat(1, intensities, intensities_build(:));

if plotfigures
    scatter(angL_build, angR_build, 'xr');
        length(directions)
    hold on
%     plot([-10 10], [-10 10])
%     plot([-10 10], [10 -10])
end


% add a few conjunctions along diagonals
angL_build = [singZeros(1), singZeros(1), singZeros(2), singZeros(2), singZeros(5), singZeros(3), singZeros(6)];
angR_build = [singZeros(1), singZeros(6), singZeros(2), singZeros(5), singZeros(2), singZeros(3), singZeros(1)];
directions_build = 999*ones(length(angL_build), 1);
intensities_build = [3, 3, 2, 2, 2, 1, 3];

angL = cat(1, angL, angL_build(:));
angR = cat(1, angR, angR_build(:));
directions = cat(1, directions, directions_build(:));
intensities = cat(1, intensities, intensities_build(:));


displacementsL = armL * tand(cat(1, 0, angL)); %this is going to be indexed 1-to-87: cannot be 90-based!!!
displacementsR = armR * tand(cat(1, 0, angR)); %this is going to be indexed 1-to-87: cannot be 90-based!!!
%moved here on Feb 5, 2019
aL87 = cat(1, 0, angL);
aR87 = cat(1, 0, angR);





% zero baseline
angL = cat(1, zeros(factor_jointZero,1), angL);
angR = cat(1, zeros(factor_jointZero,1), angR);
directions = cat(1, 999*ones(factor_jointZero,1), directions);
intensities = cat(1, zeros(factor_jointZero,1), intensities);






if plotfigures
    scatter(angL_build, angR_build);
    scatter(0, 0, '*r')
    xlabel('left antenna (deg)')
    ylabel('right antenna (deg)')
    title(sprintf('02 08 build - %d stims', length(directions)))
    grid on
    
    % % plot an iso-D curve:
    dir = 30;
    plot([cf30(dir), cf60(dir), cf120(dir)], [cf30(-dir), cf60(-dir), cf120(-dir)], '-g')
    
    % % plot a line passing for a point:
    dir = 39.5;
    plot([0, 1.7*cf60(dir)], [0, 1.7*cf60(-dir)], '-c')
    
    
    dir = 15;
    plot([cf30(dir), cf60(dir), cf120(dir)], [cf30(-dir), cf60(-dir), cf120(-dir)], '-g')
    dir = 39.5;
    plot([cf30(dir), cf60(dir), cf120(dir)], [cf30(-dir), cf60(-dir), cf120(-dir)], '-g')
    dir = 45;
    plot([cf30(dir), cf60(dir), cf120(dir)], [cf30(-dir), cf60(-dir), cf120(-dir)], '-g')
    plot([0, 1.7*cf60(dir)], [0, 1.7*cf60(-dir)], '-c')
    dir = 30;
    plot([0, 1.7*cf60(dir)], [0, 1.7*cf60(-dir)], '-c')
    dir = 15;
    plot([0, 1.7*cf60(dir)], [0, 1.7*cf60(-dir)], '-c')
    plot([-15 15], [cf60(0), cf60(0)], ':r')
    plot([-15 15], [-cf60(0), -cf60(0)], ':r')
    plot([-cf60(0), -cf60(0)], [-15 15], ':r')
    plot([cf60(0), cf60(0)], [-15 15], ':r')
    plot([cf120(0), cf120(0)], [-15 15], ':r')
    plot([-cf120(0), -cf120(0)], [-15 15], ':r')
    plot([-15 15],[-cf120(0), -cf120(0)], ':r')
    plot([-15 15],[cf120(0), cf120(0)], ':r')
    plot([-15 15],[-15 15], '--r')
    plot([-15 15],[15 -15], '--r')
    
    % % plot a single point:
    % plot(cf120(48.5), cf120(-48.5), 'x')
%     savefig('D:\Dropbox (HMS)\p2\new_2p_stimset\stimSet_20181212_f.fig')
%     export_fig('D:\Dropbox (HMS)\p2\new_2p_stimset\stimSet_20181212_wthLines_f.pdf')

    savefig('D:\Dropbox (HMS)\p2\new_2p_stimset\stimSet_20190208_f.fig')
    export_fig('D:\Dropbox (HMS)\p2\new_2p_stimset\stimSet_20190208_wthLines_f.pdf')
end

%%
seconds_addPre = 10;
seconds_addPost = 40;
window_size = 0.3*sample_fs; 
piezo90_ratio_umperVolt = 9;

% voltage_conditions = ang_displacements ./ piezo90_ratio_umperVolt;
N_conditionsTotal = length(angL) - (factor_jointZero-1); 

singleStimDuration_fsp  = sample_fs*singleStimDur;
baseline2add_start_fsp  = sample_fs*seconds_addPre;
baseline2add_end_fsp    = sample_fs*seconds_addPost;


%% fully randomize and build backbone trial-structure

% pre-build
jointZero_condition = 1; % imposed by design when building the stim vector
trialTypes = [ repmat(jointZero_condition, 1, factor_jointZero-1), 1 : N_conditionsTotal]; 
%new 1206 - same length and order as angL, unique trial label


% build backbone
both = [];
for n = 1: N_nPseudoBlocks
    both_t = repmat(trialTypes, 1, N_withinPseudoBlock);
    both = cat(2, both, both_t(randperm(length(both_t))) );
end

duration = length(both)*(singleStimDur); %seconds
duration = duration + seconds_addPre + seconds_addPost;

fprintf('------------\n')
fprintf('n stim conditions: %d\n',N_conditionsTotal)
fprintf('n repetitions: %d\n',N_withinPseudoBlock*N_nPseudoBlocks)
fprintf('single sim duration: %.2f s\n',singleStimDur)

fprintf('total duration would be: %d:%02d (min:sec) == %d sec\n', floor(duration/60), round(mod(duration, 60)), round(duration) )
% tabulate(both)
% counts = tabulate_precedingStimulus(both);


nTrials = length(both);
allOnsetsPositions = (1 : singleStimDuration_fsp : nTrials*singleStimDuration_fsp ) + baseline2add_start_fsp; % all trials, independent of type
allOffsetsPositions = (singleStimDuration_fsp : singleStimDuration_fsp : nTrials*singleStimDuration_fsp ) + baseline2add_start_fsp;


%% split backbone into left/right directions
lAnt = both; %legacy
rAnt = both; %legacy



%%%%%%%%%%%%% introduced 190207
displacementsL = displacementsL./cosd(aL87);
displacementsR = displacementsR./cosd(aR87); % the first time is to correct for big vs small angle approximation 

displacementsL = displacementsL./cosd(aL87);
displacementsR = displacementsR./cosd(aR87); % the second time is to approximate correction for piezo sensor
%%%%%%%%%%%%%


lDisplacements = displacementsL(lAnt);
rDisplacements = displacementsR(rAnt);



%% sample - build command based on DAQ proper. 
% Since it will be background acquisition, I will have to save this to a
% dat file and read and queue just a few rows per time

disp('building stim set...')
tic
lCommand = [];
rCommand = [];  % will be inverted, since this probe comes from posterior, so that stim meaning is constant

for n = 1 : length(both)

    lCommand_build = (lDisplacements(n) / piezo90_ratio_umperVolt) * ones(1, sample_fs*singleStimDur);
    lCommand = cat(2, lCommand, lCommand_build);
    
    rCommand_build = (rDisplacements(n) / piezo90_ratio_umperVolt) * ones(1, sample_fs*singleStimDur);
    rCommand = cat(2, rCommand, rCommand_build);
end

rCommand = -rCommand; %IMPORTANT

% add pre- and post- windows before filtering
rCommand = cat(2, zeros(1, baseline2add_start_fsp), rCommand, zeros(1, baseline2add_end_fsp) );
lCommand = cat(2, zeros(1, baseline2add_start_fsp), lCommand, zeros(1, baseline2add_end_fsp) );


clear lCommand_build rCommand_build
toc
%% smooth
disp('smoothing command...')
tic
commandLR(2,:) = smoothdata(rCommand, 'movmean', window_size);
commandLR(1,:) = smoothdata(lCommand, 'movmean', window_size);    
toc
% commandLR = cat(1, lCommand_f, rCommand_f);


%% write command to file
% rather than outputting the command to workspace, I need to write it to
% file.

disp('writing to file...')
identifierkey = datestr(now, 30);
logname = sprintf('p2comm_%s.bin', identifierkey); % will save in current folder wich should be runfolder

fid = fopen(logname, 'w');
fwrite(fid, commandLR, 'single');
fclose(fid);



% % note that I can read it back, specifying 2 rows and as many columns I
% % want my chunk to have. No need to specify the file position, which
% % continues from where the previous read ended, but only if file is never
% % closed and reopened. Example:
% fid = fopen(logname, 'r+');
% A = fread(fid, [2, chuncklengt] ,'single' );
% % when it reaches the end, A is empty.


% 
% %% let's make a simpler example first
% saveData = [1:20; 101:120];
% fid = fopen(logname, 'w');
% fwrite(fid, saveData, 'single');
% fclose(fid);
% 
% 
% %% try to read back chunks
% fid = fopen(logname, 'r+');
% A = fread(fid, [2,5], 'single') %works indeed
% 
% % repeat to see if continues automatically from previous position?
% A = fread(fid, [2,5], 'single') % it does indeed
% 
% % check if this syntax is equivalent
% fseek(fid,0,'cof');
% A = fread(fid, [2,5], 'single') % it is indeed

% fclose(fid)
%% metadata to workspace and saved to disk
disp('preparing and saving metadata output...')
stimset.key.ID = identifierkey;
stimset.key.logname = logname;

stimset.joint.trialIndices = both;
stimset.joint.allOnsetsPositions = allOnsetsPositions;
stimset.joint.allOffsetsPositions = allOffsetsPositions;

stimset.singleAnt.L_indices = lAnt;
stimset.singleAnt.R_indices = rAnt;
stimset.singleAnt.L_displacements = lDisplacements; %positive is always push toward head, negative pull
stimset.singleAnt.R_displacements = rDisplacements; %positive is always push toward head, negative pull

%stimset.dec contains redundant information, that may simplify things
%later, possibly
stimset.dec.jointZero_index = jointZero_condition;
stimset.dec.trialTypes = trialTypes;
stimset.dec.angL = angL;
stimset.dec.angR = angR;
stimset.dec.armL = armL;
stimset.dec.armR = armR;
stimset.dec.directions = directions;
stimset.dec.intensities = intensities;

stimset.dec.i2displ_model = samples(6,:);
stimset.dec.i2displ_decoder = iWind; 


stimset.userinput.N_withinPseudoBlock   = N_withinPseudoBlock;
stimset.userinput.factor_jointZero      = factor_jointZero;
stimset.userinput.N_nPseudoBlocks       = N_nPseudoBlocks;
stimset.userinput.singleStimDur         = singleStimDur;
stimset.userinput.DAQ_fs                = sample_fs;

stimset.experimentHandle    = experimentHandle;   % for completeness of info all in the same place, and because this info is not included in the metadata filename


% save to disk
save(sprintf('metadata_%s.mat', identifierkey), '-struct', 'stimset' );

disp('Done.')

end

