load('/Users/galileo/Dropbox (HMS)/p2/dataTable_p2.mat', 'T')


t = 1;

flyNum     = T.flynum(t);
cellNum    = 1;
disp(T.metadatafiles{t})

%select which datafiles
datafiles = uipickfiles('FilterSpec', sprintf('/Users/galileo/Dropbox (HMS)/p2/fly%3d_PP/data_*', flyNum));
disp('loading...')
clear data
for i = 1:length(datafiles)
    disp(i)
    data(i) = load(datafiles{i});
end
disp('done.')

aL = data(1).aL;
aR = data(1).aR;

%% pull out single stim example traces
figFolder = fullfile(fileparts(T.linkToDataFile{t}), sprintf('singleStimsTraces_position_%s', data(i).metadata.key.ID));
mkdir(figFolder);


i = 1;

for stim = 1:length(data(i).perist)
    
    stim_I = find(data(i).trialIndicesFullBlock == stim);
    ps = -2;
    clear prePost_stimPositions
    for p = 1:5
        prePost_stimPositions(p,:) = stim_I + ps;
        ps = ps+1;
    end
    prePost_stimPositions(prePost_stimPositions<1) = nan;
    prePost_stimPositions(prePost_stimPositions>length(data(i).trialIndicesFullBlock)) = nan;
    use_PPsP = prePost_stimPositions(:);
    
    prePost_stimArray = nan(size(prePost_stimPositions));
    prePost_stimArray(~isnan(use_PPsP)) = data(i).trialIndicesFullBlock(use_PPsP(~isnan(use_PPsP)));
    prePost_stimArray = reshape(prePost_stimArray,size(prePost_stimPositions));
    
    prePost_stimArray = prePost_stimArray';
    
    figure('Color', [1 1 1]) ; hold on
    for trial = 1:length(stim_I)
        for s = 1:5
            L(s) = aL(st2rep)
            R(s) = aR(st2rep)
            st2rep = prePost_stimArray(trial,s);
            subplot(length(stim_I), 5, (trial - 1)*5 + s)
            if ~isnan(st2rep)
                scatter(aL(st2rep), aR(st2rep), 'k', 'filled')
            end
            xlim([-15, 15])
            ylim([-15, 15])
            axis square
            box on
            ax = gca;
            ax.YTick = [];
            ax.XTick = [];
            pause
        end
    end
    stringName = sprintf('%s_%s_position_stim%02d', data(i).metadata.key.ID(1:8), data(i).metadata.key.ID(10:end), stim);
    export_fig(fullfile(figFolder, sprintf('%s.eps',stringName)))
    close

end

%%
figure; hold on; 
for i = 1:5
    scatter(aL(prePost_stimArray(9,i)), aR(prePost_stimArray(9,i)), 'filled')
end