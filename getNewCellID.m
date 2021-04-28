function experimentHandle = getNewCellID(experimentHandle)

%% ask me what to do, and do it
lastcell_basename = sprintf('fly%03d_cell%02d',experimentHandle.flyNum, experimentHandle.cellNum);
answer = questdlg(sprintf('Last cell was: %s. \nWhere to?', lastcell_basename), ...
	'new acquisition location', ...
	'continue with last cell','record from a new cell','continue with last cell');
% Handle response
if strcmp(answer, 'record from a new cell')
        experimentHandle.cellNum = experimentHandle.cellNum +1;
        experimentHandle.basename = sprintf('fly%03d_cell%02d',experimentHandle.flyNum, experimentHandle.cellNum);
        experimentHandle.cellfolder = fullfile(experimentHandle.flyfolder, experimentHandle.basename);
        mkdir(experimentHandle.cellfolder);
        experimentHandle.recycleMetadata = 1;
end

disp(experimentHandle)
cd(experimentHandle.cellfolder)

end