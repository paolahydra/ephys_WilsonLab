function experimentHandle = getNewFlyID(rootfolder, lastFly)
%experimentHandle = getNewFlyID(rootfolder);
% lets you choose between last cell, new cell or new fly
% assumes that folder naming  for fly- and cell- folders in rootfolder is fixed and reserved
% required input: rootfolder
% optional input: lastFly. This will override search for the last flyfolder
% created and will allow you to make a new cell in the specified flyfolder
% created and point 
%
% output fields:
% 'flyNum'   
% 'cellNum'   
% 'basename' 
% 'flyfolder'
% 'cellfolder'
%
% 
% PP 03/15/2018

narginchk(1,2)

%% find actual lastFlyNum
expstr =  'fly\d{3}_PP';
recursive = 0;
returnfullpath = 0;
dironly = 1;
dirlist = regexpdir(rootfolder, expstr, recursive, returnfullpath, dironly);

for i = 1:length(dirlist)
    currentflies(i,1) = str2num(dirlist{i}(4:6));
end
currentflies = sort(currentflies);
lastFlyNum = currentflies(end);


if ~exist('lastFly','var')
    lastFly = lastFlyNum;
end



%% get all current cells of last fly
lastFlyFolder = fullfile(rootfolder, sprintf('fly%03d_PP',lastFly));

expstr =  'fly\d{3}_cell\d{2}';
recursive = 0;
returnfullpath = 0;
dironly = 1;
celldirlist = regexpdir(lastFlyFolder, expstr, recursive, returnfullpath, dironly);

currentcells = [];
for i = 1:length(celldirlist)
    currentcells(i,1) = str2num(celldirlist{i}(end-1:end));
end
if ~isempty(currentcells)
    currentcells = sort(currentcells);
    lastCell_lastFly = currentcells(end);
else
    lastCell_lastFly = 0;
end

%% ask me what to do, and do it
lastcell_basename = sprintf('fly%03d_cell%02d',lastFly, lastCell_lastFly);
answer = questdlg(sprintf('Last cell was: %s. \nWhere to?', lastcell_basename), ...
	'new acquisition location', ...
	'last cell','new cell','new fly','new fly');
% Handle response
switch answer
    case 'last cell'
        experimentHandle.flyNum = lastFly;
        experimentHandle.cellNum = lastCell_lastFly;
        experimentHandle.basename = sprintf('fly%03d_cell%02d',lastFly, lastCell_lastFly);
        experimentHandle.flyfolder = lastFlyFolder;
        experimentHandle.cellfolder = fullfile(lastFlyFolder, experimentHandle.basename);
        assert(lastCell_lastFly~=0, 'Error: lastCell folder does not exist.')
    case 'new cell'
        experimentHandle.flyNum = lastFly;
        experimentHandle.cellNum = lastCell_lastFly +1;
        experimentHandle.basename = sprintf('fly%03d_cell%02d',lastFly, experimentHandle.cellNum);
        experimentHandle.flyfolder = lastFlyFolder;
        experimentHandle.cellfolder = fullfile(lastFlyFolder, experimentHandle.basename);
        mkdir(experimentHandle.cellfolder);
        
    case 'new fly'
        experimentHandle.flyNum = lastFlyNum +1; %note this will increase no matter the overriding input
        experimentHandle.cellNum = 1;
        experimentHandle.basename = sprintf('fly%03d_cell%02d',experimentHandle.flyNum, experimentHandle.cellNum);
        experimentHandle.flyfolder = fullfile(rootfolder, sprintf('fly%03d_PP', experimentHandle.flyNum) );
        experimentHandle.cellfolder = fullfile(experimentHandle.flyfolder, experimentHandle.basename);
        mkdir(experimentHandle.cellfolder);

end
experimentHandle.recycleMetadata = 0;
disp(experimentHandle)
cd(experimentHandle.cellfolder)
end