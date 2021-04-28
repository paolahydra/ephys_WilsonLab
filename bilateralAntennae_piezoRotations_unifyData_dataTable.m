load('/Users/galileo/Dropbox (HMS)/p2/dataTable_p2.mat', 'T')
% % this is the main code, but still some blocks are only left in  
% edit bilateralAntennae_piezoRotations_unifyData_fromSingleRuns_May

singleAntData = matfile('/Users/galileo/Dropbox (HMS)/p2/singleAntennaData.mat', 'Writable', true);
Nreps = 1e5;
sampleZeroMean(Nreps) = nan;

%% load and resave dataTable

for t = 1
    
    flyNum     = T.flynum(t);
    if flyNum-1000>0
        cellNum = mod(flyNum, 10);
        flyNum = floor(flyNum/10);
    else
        cellNum    = 1;
    end
    
    disp(T.metadatafiles{t})
    
    datafiles = uipickfiles('FilterSpec', sprintf('/Users/galileo/Dropbox (HMS)/p2/fly%3d_PP/data_fly%3d_cell%02d_*', flyNum, flyNum, cellNum));
    disp('loading...')
    clear data
    for i = 1:length(datafiles)
        disp(i)
        data(i) = load(datafiles{i});
    end
    disp('done.')
    singleAntData.datafiles(1, t) = {datafiles};
    
%     meanFR = cat(1, data.meanFR);             % this was wrong, i.e. not weighted by block size
%     meanfiltVm = cat(1, data.meanfiltVm);
%     meanFR = mean(meanFR,1);
%     meanfiltVm = mean(meanfiltVm,1);
    meanFR = [];
    meanfiltVm = [];
    for st = 1:length(data(1).meanFR) % I am tired. I'll do it the safe way
        a = [];
        v = [];
        for i = 1:length(datafiles)
            a = cat(1, a, data(i).FR{st});
            v = cat(1, v, data(i).filtVm{st});
        end
        meanFR(st) = mean(a);
        meanfiltVm(st) = mean(v);
    end
    
    aL = data(1).aL;
    aR = data(1).aR;
    intensities = data(1).intensities;
    directions = data(1).directions;
    
    linkToDataFile = fullfile(fileparts(datafiles{1}), sprintf('dataFile_fly%03d_cell%02d.mat',flyNum, cellNum));
    save(linkToDataFile, 'meanFR', 'meanfiltVm', 'aL', 'aR', 'intensities', 'directions');
    
    % do not build dataTable up (too big!)
    T.linkToDataFile{t} = linkToDataFile;
    
    if ~ismember(t, 8:10)
        I30 = intensities == 30; %patch
        I60 = intensities == 60;
        I120 = intensities == 120;
        
        T.meanFR_I30(t,:) = meanFR(I30)-meanFR(1);
        T.meanFR_I60(t,:) = meanFR(I60)-meanFR(1);
        T.meanFR_I120(t,:) = meanFR(I120)-meanFR(1);
        
        T.meanVm_I30(t,:) = meanfiltVm(I30)-meanfiltVm(1);
        T.meanVm_I60(t,:) = meanfiltVm(I60)-meanfiltVm(1);
        T.meanVm_I120(t,:) = meanfiltVm(I120)-meanfiltVm(1);
        
        T.directions_I30(t,:) = directions(I30);
        T.directions_I60(t,:) = directions(I60);
        T.directions_I120(t,:) = directions(I120);
    end
    clc
    save('/Users/galileo/Dropbox (HMS)/p2/dataTable_p2.mat', 'T', '-append');
    % trialIndicesFullBlock = data(1).trialIndicesFullBlock;
    
    
    scatfold = '/Users/galileo/Dropbox (HMS)/p2/correctedScatterplots_190603';
    
    %% dir tuning
%     if ~ismember(t, 8:10)
%         f = figure; hold on
%         plot(directions(I30), meanFR(I30)-meanFR(1));
%         plot(directions(I60), meanFR(I60)-meanFR(1));
%         plot(directions(I120), meanFR(I120)-meanFR(1));
%         legend('30', '60', '120')
%         legend boxoff
%         legend('Location', 'NorthWest')
%         xlabel('estimated wind direction')
%         ylabel('mean firing rate change with respect to rest')
%         title(T.flynum(t), 'Interpreter', 'none')
%         
%         savefig(fullfile(scatfold, sprintf('AVG_dirTuning_FR_RestSubtr_%d', T.flynum(t))))
%         export_fig(fullfile(scatfold, sprintf('AVG_dirTuning_FR_RestSubtr_%d.pdf', T.flynum(t))))
%     end
    
    %% general color-coded scatter
    
    figure; hold on; axis image
    xlabel(' Left antenna (ipsi) - angular displ')
    ylabel(' Right antenna (contra) - angular displ')
    scatter(aL, aR, 200, meanfiltVm-meanfiltVm(data(1).metadata.dec.jointZero_index), 'filled', 'MarkerEdgeColor', 'k')
    title(sprintf('Vm (N = %d)\n fly %d', T.NtrialsIncluded(t), T.flynum(t)), 'Interpreter', 'none')
    colormap(bluewhitered(256)), cb = colorbar;
    cb.TickLabels = num2str((cb.Ticks+meanfiltVm(data(1).metadata.dec.jointZero_index))', '%2.1f');
    cb.Label.String = '(mV)';
    xlim([-15 15])
    ylim([-15 15])
    set(gca, 'TickDir', 'out')
    export_fig(fullfile(scatfold, sprintf('dtAVG_Vm_scatter_fly%d.pdf', T.flynum(t))))
    
    
    figure; hold on; axis image
    xlabel(' Left antenna (ipsi) - angular displ')
    ylabel(' Right antenna (contra) - angular displ')
    scatter(aL, aR, 200, meanFR-meanFR(data(1).metadata.dec.jointZero_index), 'filled', 'MarkerEdgeColor', 'k')
    title(sprintf('firing rate (N = %d)\n fly %d', T.NtrialsIncluded(t), T.flynum(t)), 'Interpreter', 'none')
    colormap(bluewhitered(256)), cb = colorbar;
    cb.TickLabels = num2str((cb.Ticks+meanFR(data(1).metadata.dec.jointZero_index))', '%2.1f');
    cb.Label.String = '(Hz)';
    xlim([-15 15])
    ylim([-15 15])
    set(gca, 'TickDir', 'out')
    export_fig(fullfile(scatfold, sprintf('dtAVG_FR_scatter_fly%d.pdf', T.flynum(t))))
    
    
    %% interpolate scatter FR (natural)
    
    % bring aL and aR in mash grid format
    V = meanFR'-meanFR(data(1).metadata.dec.jointZero_index);
    F = scatteredInterpolant(aL(:),aR(:),V);
    % F.Method = 'linear';
    % ExtrapolationMethod = 'nearest';
    F.Method = 'natural';
    ExtrapolationMethod = 'nearest';
    
    [Xq,Yq] = meshgrid(min(aL):0.25:max(aL));
    Vq = F(Xq, Yq);
    
    figure;
    xlabel('L (deg)')
    ylabel('R (deg)')
    h = surf(Xq,Yq,Vq);
    
    az = 0;
    el = 90;
    view(az, el);
    
    % ylim([min(aL), max(aL)])
    % xlim([min(aL), max(aL)])
    xlim([-15 15])
    ylim([-15 15])
    axis square
    a = gca;
    a.XAxis.Visible = 'off';
    a.YAxis.Visible = 'off';
    a.XGrid = 'off';
    a.YGrid = 'off';
    h.EdgeColor = 'none';
    
    colormap(bluewhitered(256)), cb = colorbar;
    cb.TickLabels = num2str((cb.Ticks+meanFR(data(1).metadata.dec.jointZero_index))', '%2.1f');
    cb.Label.String = '(Hz)';
    
    title(data(1).metadata.experimentHandle.basename, 'Interpreter', 'none')
    
    % export_fig(sprintf('dtAVG_interpSc_FR_natInt_NearestExt_fly%d.jpg', T.flynum(t)))
    export_fig(fullfile(scatfold, sprintf('dtAVG_interpSc_FR_natInt_NearestExt_fly%d.pdf', T.flynum(t))))
    close


    %% single antenna contribution and CI - saving data - 190603
    if ~ismember(t, 8:10)
        ixL = find(aR==0);
        ixL = ixL([2:4,1,5:7]);
        aL_ixL = aL(ixL);
        
        ixR = find(aL==0);
        ixR = ixR([2:4,1,5:7]);
        aR_ixR = aR(ixR);  

        FRZero = [];
        FRStL = [];
        FRStR = [];
        for bl = 1:length(data)
            %     assert(data(bl).metadata.userinput.N_withinPseudoBlock == 1, 'more than one rep per block not exploited')
            FRZero = cat(1, FRZero, data(bl).FR{1});
            FRStL = cat(1, FRStL, cat(2,data(bl).FR{setdiff(ixL,1)}) );
            FRStR = cat(1, FRStR, cat(2,data(bl).FR{setdiff(ixR,1)}) );
        end
        if length(directions) == 84
            FRZero = reshape(FRZero, 7, []);
        else
            FRZero = mean(reshape(FRZero, 4, []));
        end
        sz1 = size(FRZero, 1);
        sz2 = size(FRZero, 2);
        
        singleAntData.FRZero(1,t) = {FRZero};
        singleAntData.FRStL(1,t) = {FRStL};
        singleAntData.FRStR(1,t) = {FRStR};
        
        %% confidence interval based on SE:
        singleAntData.CI95(1,t) = {[mean(FRZero(:)) - std(mean(FRZero))/sqrt(sz2) * 1.96, ...
                    mean(FRZero(:)) + std(mean(FRZero))/sqrt(sz2) * 1.96]};
                
        %% shuffle within columns of FRZero to get random samples of temporally
        %sorted samples
        
        FRZero = FRZero(:); %still 7 elements per block
        addendum = 0:sz1:numel(FRZero)-1;
        
        tic
        for i = 1:Nreps
            a = randi(sz1, 1,sz2);
            idx = a + addendum;
            sampleZeroMean(i) = mean(FRZero(idx));
        end
        toc
        singleAntData.CI95bootstrap(1,t) = {[prctile(sampleZeroMean, 2.5),...
            prctile(sampleZeroMean,97.5)]};
        
        %% show it
        figure; hold on
        FR_L_only = meanFR(ixL);
        FR_R_only = meanFR(ixR);
        h1 = plot(aL_ixL, FR_L_only', '-<');
        h2 = plot(aR_ixR, FR_R_only', '->');
        
        zero = aL_ixL==0;
        CI_centered = singleAntData.CI95bootstrap(1,t);
        CI_centered = CI_centered{1} - FR_L_only(zero);
        plot([aL_ixL(1), aL_ixL(end)], [meanFR(ixL(4)), meanFR(ixL(4))], '-k')
        plot([aL_ixL(1), aL_ixL(end)], [meanFR(ixL(4))+CI_centered(1), meanFR(ixL(4))+CI_centered(1)], ':k')
        plot([aL_ixL(1), aL_ixL(end)], [meanFR(ixL(4))+CI_centered(2), meanFR(ixL(4))+CI_centered(2)], ':k')
        legend([h1, h2],{'L moving, R at rest', 'R moving, L at rest'}, 'Location', 'best')
        title(flyNum)
        export_fig(fullfile(scatfold, sprintf('AVG_singleAntennaeContribution_fly%d.pdf', flyNum)))

        
    else
        aL_ixL = aL(aR==0); %these are already sorted!
        aR_ixR = aR(aL==0);  %these are already sorted!
        
        ixL = find(aR==0);
        ixR = find(aL==0);
        
        FRZero = [];
        FRStL = [];
        FRStR = [];
        for bl = 1:length(data)
            %     assert(data(bl).metadata.userinput.N_withinPseudoBlock == 1, 'more than one rep per block not exploited')
            FRZero = cat(1, FRZero, data(bl).FR{data(1).metadata.dec.jointZero_index});
            FRStL = cat(1, FRStL, cat(2,data(bl).FR{setdiff(ixL,data(1).metadata.dec.jointZero_index)}) );
            FRStR = cat(1, FRStR, cat(2,data(bl).FR{setdiff(ixR,data(1).metadata.dec.jointZero_index)}) );
        end

        FRZero = reshape(FRZero, data(1).metadata.userinput.factor_jointZero*data(1).metadata.userinput.N_withinPseudoBlock, []);
        
        sz1 = size(FRZero, 1);
        sz2 = size(FRZero, 2);
        
        singleAntData.FRZero(1,t) = {FRZero};
        singleAntData.FRStL(1,t) = {FRStL};
        singleAntData.FRStR(1,t) = {FRStR};
        
        %% confidence interval based on SE: % this is ot fear because sampling size of zero is always larger
        singleAntData.CI95(1,t) = {[mean(FRZero(:)) - std(mean(FRZero))/sqrt(sz2) * 1.96, ...
                    mean(FRZero(:)) + std(mean(FRZero))/sqrt(sz2) * 1.96]};
                
        %% shuffle within columns of FRZero to get random samples of temporally
        %sorted samples
        FRZero = FRZero(:); %still sz1 elements per block
        addendum = 0:sz1:numel(FRZero)-1;
        
        tic
        for i = 1:Nreps
            a = randi(sz1, 1,sz2);
            idx = a + addendum;
            sampleZeroMean(i) = mean(FRZero(idx));
        end
        toc
        singleAntData.CI95bootstrap(1,t) = {[prctile(sampleZeroMean, 2.5),...
            prctile(sampleZeroMean,97.5)]};
        
        %% show it
        figure; hold on
        FR_L_only = meanFR(ixL);
        FR_R_only = meanFR(ixR);
        h1 = plot(aL_ixL, FR_L_only', '-<');
        h2 = plot(aR_ixR, FR_R_only', '->');
        
        zero = aL_ixL==0;
        CI_centered = singleAntData.CI95bootstrap(1,t);
        CI_centered = CI_centered{1} - FR_L_only(zero);
        plot([aL_ixL(1), aL_ixL(end)], [meanFR(ixL(4)), meanFR(ixL(4))], '-k')
        plot([aL_ixL(1), aL_ixL(end)], [meanFR(ixL(4))+CI_centered(1), meanFR(ixL(4))+CI_centered(1)], ':k')
        plot([aL_ixL(1), aL_ixL(end)], [meanFR(ixL(4))+CI_centered(2), meanFR(ixL(4))+CI_centered(2)], ':k')
        legend([h1, h2],{'L moving, R at rest', 'R moving, L at rest'}, 'Location', 'best')
        title(flyNum)
        export_fig(fullfile(scatfold, sprintf('AVG_singleAntennaeContribution_fly%d_cell%d.pdf', flyNum, cellNum)))


    end

end
save('/Users/galileo/Dropbox (HMS)/p2/dataTable_p2.mat', 'T', '-append');