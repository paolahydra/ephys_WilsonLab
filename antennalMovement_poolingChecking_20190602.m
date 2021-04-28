load('/Users/galileo/Dropbox (HMS)/p2/dataTable_p2.mat', 'T')

for t = 1:9
    clc
    flyNum     = T.flynum(t);
    fly(t).AI = load(sprintf('/Users/galileo/Dropbox (HMS)/p2/antImages_fly%d.mat',flyNum));
end

% clc
% for t = 1:7
%     fprintf('fly %d\n\tarmL:\t%.1f\n\tarmR:\t%.1f\n', T.flynum(t), fly(t).AI.calibrationData(end).armL, fly(t).AI.calibrationData(end).armR);
% end

% fly 292
% 	armL:	237.8
% 	armR:	261.4
% fly 290
% 	armL:	195.0
% 	armR:	235.3
% fly 267
% 	armL:	187.5
% 	armR:	253.2
% fly 264
% 	armL:	175.0
% 	armR:	196.0
% fly 263
% 	armL:	231.4
% 	armR:	244.4
% fly 279
% 	armL:	146.0
% 	armR:	177.0
% fly 283
% 	armL:	156.3
% 	armR:	150.5

%%
fly2reg = 9;
cam = 2;
i = 5; %ordinal mumber pof piezoDispl file
fr = 9;
registrandum = fly(fly2reg).AI.A{i}(cam).frames(:,:,fr);
figure; imshow(registrandum);

%%
figure; imshow(fly(8).AI.A{i}(cam).frames(:,:,fr));

%% 
clear refFrame refFrame_pushpull
for t = 1:8
    refFrame(:,:,t) = fly(t).AI.A{end}(cam).frames(:,:,fr);
    pushpull = imshowpair(fly(t).AI.A{end}(cam).frames(:,:,4), fly(t).AI.A{end}(cam).frames(:,:,5)); %smaller stims are 24 um
    refFrame_pushpull(:,:,:,t) = pushpull.CData;
end
refFrame(:,:,fly2reg) = registrandum;
registrandum_pushpull = imshowpair(fly(fly2reg).AI.A{i}(cam).frames(:,:,7), fly(fly2reg).AI.A{i}(cam).frames(:,:,8)); %largest stims are 20um
refFrame_pushpull(:,:,:,fly2reg) = registrandum_pushpull.CData;
%% show all
% for t = 1:fly2reg %show arm
%     figure(t)
%     imshow(refFrame(:,:,t));
%     title(sprintf('%d', T.flynum(t)), 'Interpreter', 'none')
% end
for t = 1:fly2reg 
    figure(t)
    imshow(refFrame_pushpull(:,:,:,t));
    title(sprintf('%d', T.flynum(t)), 'Interpreter', 'none')
end