load('/Users/galileo/Dropbox (HMS)/p2/dataTable_p2.mat', 'T')


for t = 9
    clc
    flyNum     = T.flynum(t);
    if 1000-flyNum>0
        cellNum    = 1;
        flyNum3digits = flyNum;
    else
        cellNum = mod(flyNum, 10);
        flyNum3digits = floor(flyNum/10);
    end
        
    disp(T.metadatafiles{t})
    metadataFiles = dir(sprintf('/Users/galileo/Dropbox (HMS)/p2/fly%3d_PP/fly%3d_cell%02d/metadata_*.mat', flyNum3digits, flyNum3digits, cellNum));
    cellFolder = metadataFiles(1).folder;
	metadataFiles = {metadataFiles(:).name};
    antRecFiles = dir(sprintf('/Users/galileo/Dropbox (HMS)/p2/fly%3d_PP/piezosDispl_fly%3d*.mat', flyNum3digits, flyNum3digits));
    flyFolder = antRecFiles(1).folder;
    antRecFiles = {antRecFiles(:).name};
    
    chronoFiles = cat(2, metadataFiles, antRecFiles); %not as important!!
    for c = 1:length(chronoFiles)
        cx = strfind(chronoFiles{c}, '_201');
        numFiles{c} = chronoFiles{c}(cx+1:end);
    end
    [~, I] = sort(numFiles);
    clear numFiles c cx
    disp('chronologically sorted files:')
    chronoFiles = chronoFiles(I);
    disp(chronoFiles')
    
    calibrationI = false(size(antRecFiles));
    calibrationData = [];
    clear A
    for i = 1:length(antRecFiles)
        a = load(fullfile(flyFolder, antRecFiles{i}));
        if isfield(a, 'dispFrames')
            calibrationI(i) = true;
            calibrationData(i).savedCalData = a;
            calibrationData(i).armL = nanmedian(calibrationData(i).savedCalData.armL_um);
            calibrationData(i).armR = nanmedian(calibrationData(i).savedCalData.armR_um);
        end
        A{i} = a.allFr;
    end  
    save(sprintf('/Users/galileo/Dropbox (HMS)/p2/antImages_fly%d.mat', flyNum), 'A', 'calibrationI', 'calibrationData', ...
    'chronoFiles', 'antRecFiles')
end

% save(sprintf('/Users/galileo/Dropbox (HMS)/p2/antImages_fly%3d.mat', flyNum), 'A', 'calibrationI', 'calibrationData', ...
%     'chronoFiles', 'antRecFiles')


%% %% mostra uno stesso frame over time (across all recs)

offsetFig = 10;
warning('off')

cam = 2;    % 2 is bottom camera
fr = 9;     % 3 6 9 are rest, piezo-on
% cab(1)
for i = 1:length(A)
    figure(i+offsetFig)
    imshow(A{i}(cam).frames(:,:,fr));
    title(sprintf('frame %d -- rec %d: %s',fr, i, antRecFiles{i}), 'Interpreter', 'none')
end

%% %% mostra tutti i frames di uno stesso recording

offsetFig = 10;
warning('off')

cam = 3;    % 2 is bottom camera
% length(A)
i = 4;     
% cab(1)
for fr = 1:9 % 3 6 9 are rest, piezo-on - 7 is biggest PULL
    figure(fr+offsetFig)
    imshow(A{i}(cam).frames(:,:,fr));
    title(sprintf('frame %d -- rec %d: %s',fr, i, antRecFiles{i}), 'Interpreter', 'none')
end

