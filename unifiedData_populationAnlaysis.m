%post-unification script!
scatfold = '/Users/galileo/Dropbox (HMS)/p2/correctedScatterplots_190603';
    
load('/Users/galileo/Dropbox (HMS)/p2/dataTable_p2.mat', 'T')
singleAntData = matfile('/Users/galileo/Dropbox (HMS)/p2/singleAntennaData.mat');
% load stuff up. Assume X vars are different
for t = 1:size(T,1)
    
    mfile = matfile(T.linkToDataFile{t});
    
    data(t).flyNum = T.flynum(t);
    data(t).aL = mfile.aL;
    data(t).aR = mfile.aR;
    data(t).dir = mfile.directions;
    data(t).int = mfile.intensities;
    data(t).Vm = mfile.meanfiltVm;
    data(t).FR = mfile.meanFR;
end
clear mfile
% I30 = intensities == 30; 
% I60 = intensities == 60;
% I120 = intensities == 120;

%% population direction tuning curve at highest intensity
for t = 1:7
    aR = data(t).aR;
    aL = data(t).aL;
    directions = data(t).dir;
    intensities = data(t).int;
    meanFR = data(t).FR;
    
    I30 = intensities == 30; %patch
    I60 = intensities == 60;
    I120 = intensities == 120;
    
    unifSingAntData(t).dirI30 = meanFR(I30)-meanFR(1);
    unifSingAntData(t).dirI60 = meanFR(I60)-meanFR(1);
    unifSingAntData(t).dirI120 = meanFR(I120)-meanFR(1);
end    

figure; hold on
h = plot(directions(I120), cat(1,unifSingAntData.dirI120), '-o', 'MarkerSize', 4);
set(h, {'MarkerFaceColor'}, get(h,'Color')); 
legend(num2str(T.flynum))
legend boxoff
legend('Location', 'SouthEast')
ylabel('FR (Hz)')

savefig(fullfile(scatfold, sprintf('popul_dirTuning_FR_RestSubtr_I120')))
export_fig(fullfile(scatfold, sprintf('popul_dirTuning_FR_RestSubtr_I120.pdf')))


% normalize between -1 and 1 first
Y = cat(1,unifSingAntData.dirI120);
Y = bsxfun(@rdivide, Y, max(Y,[],2) );
NegPeakDivide = max(ones(size(Y,1), 1), abs(min(Y,[],2)) );
Y = bsxfun(@rdivide, Y, NegPeakDivide );

figure; hold on
h = plot(directions(I120), Y, '-o', 'MarkerSize', 4);
set(h, {'MarkerFaceColor'}, get(h,'Color')); 
legend(num2str(T.flynum))
legend boxoff
legend('Location', 'SouthEast')
ylabel('norm -1 to 1 FR')

savefig(fullfile(scatfold, sprintf('popul_dirTuning_normFR_I120')))
export_fig(fullfile(scatfold, sprintf('popul_dirTuning_normFR_I120.pdf')))


% smooth a bit?
Y = cat(1,unifSingAntData.dirI120);
Y = movmean(Y, 3, 2);
Y = bsxfun(@rdivide, Y, max(Y,[],2) );
NegPeakDivide = max(ones(size(Y,1), 1), abs(min(Y,[],2)) );
Y = bsxfun(@rdivide, Y, NegPeakDivide );

figure; hold on
h = plot(directions(I120), Y, '-o', 'MarkerSize', 4);
set(h, {'MarkerFaceColor'}, get(h,'Color')); 
legend(num2str(T.flynum))
legend boxoff
legend('Location', 'SouthEast')
ylabel('norm -1 to 1: smoothed(FR), movmean 3')

savefig(fullfile(scatfold, sprintf('popul_dirTuning_normsmoothFR_I120')))
export_fig(fullfile(scatfold, sprintf('popul_dirTuning_normsmoothFR_I120.pdf')))

%% sort out single antennae contribution
for t = 1:size(T,1)
    aR = data(t).aR;
    aL = data(t).aL;
    directions = data(t).dir;
    intensities = data(t).int;
    meanFR = data(t).FR;   
    
    if ~ismember(t, 8:10)
        ixL = find(aR==0);
        ixL = ixL([2:4,1,5:7]);
        aL_ixL = aL(ixL);
        
        ixR = find(aL==0);
        ixR = ixR([2:4,1,5:7]);
        aR_ixR = aR(ixR);   % xaxis - intervals are not linear, plus values might not be sorted? - but they are
        
        unifSingAntData(t).aL_only = aL_ixL;
        unifSingAntData(t).aR_only = aR_ixR;
        unifSingAntData(t).FR_L_only = meanFR(ixL);
        unifSingAntData(t).FR_R_only = meanFR(ixR);
        
    else % show single antenna contribution for 223 233 flies      
        aL_ixL = aL(aR==0); %these are already sorted!
        aR_ixR = aR(aL==0);  %these are already sorted!
        
        unifSingAntData(t).aL_only = aL_ixL;
        unifSingAntData(t).aR_only = aR_ixR;
        unifSingAntData(t).FR_L_only = meanFR(aR==0);
        unifSingAntData(t).FR_R_only = meanFR(aL==0);
        
    end
end
clearvars -except unifSingAntData data T singleAntData scatfold

%% plot all lines for each antenna
% figure; hold on
% for t = 1:size(T,1)
%     plot(unifData(t).aL_only, unifData(t).FR_L_only);
% end                                                 % all over the place

figure; hold on
a = gca;
a.ColorOrder = cat(1, a.ColorOrder, [38 34 97; 235 0 139; 0 104 56]./255);
for t = 1:size(T,1)
    zero = unifSingAntData(t).aL_only==0;
    plot(unifSingAntData(t).aL_only, unifSingAntData(t).FR_L_only - unifSingAntData(t).FR_L_only(zero));
end
title('L antenna')
legend(num2str(T.flynum))
legend('Location', 'Best')
ylim([-8 10])
export_fig(fullfile(scatfold, sprintf('allFlies_singleAntenna_LEFTmoving.pdf')))


figure; hold on
a = gca;
a.ColorOrder = cat(1, a.ColorOrder, [38 34 97; 235 0 139; 0 104 56]./255);
for t = 1:size(T,1)
    zero = unifSingAntData(t).aL_only==0; %L-R does not matter, there is just one zero.
    plot(unifSingAntData(t).aR_only, unifSingAntData(t).FR_R_only - unifSingAntData(t).FR_R_only(zero));
end
title('R antenna')
legend(num2str(T.flynum))
legend('Location', 'Best')
ylim([-8 10])
export_fig(fullfile(scatfold, sprintf('allFlies_singleAntenna_RIGHTmoving.pdf')))

%% plot one cell at the time, with CI of zero - some error, but already done in bilateralAnt_unify_datatable...
% for t = 9:size(T,1)
%     figure; hold on
%     title(T.flynum(t))
%     zero = unifData(t).aL_only==0;
%     CI_centered = singleAntData.CI95bootstrap(1,t);
%     CI_centered = CI_centered{1} - unifData(t).FR_L_only(zero);
%     plot(unifData(t).aL_only, unifData(t).FR_L_only - unifData(t).FR_L_only(zero), '-<');
%     plot(unifData(t).aL_only, unifData(t).FR_R_only - unifData(t).FR_R_only(zero),  '->');
%     plot([unifData(t).aL_only(1), unifData(t).aL_only(end)], [unifData(t).FR_L_only(zero)+CI_centered(2), unifData(t).FR_L_only(zero)+CI_centered(2)], ':k')
%     plot([unifData(t).aL_only(1), unifData(t).aL_only(end)], [unifData(t).FR_L_only(zero)+CI_centered(1), unifData(t).FR_L_only(zero)+CI_centered(1)], ':k')
% end


%% averaging push and averaging pull displacements
for t = 1:size(T,1)
    zeroI = find(unifSingAntData(t).aL_only==0);
    FR_L = unifSingAntData(t).FR_L_only - unifSingAntData(t).FR_L_only(zeroI);
    FR_R = unifSingAntData(t).FR_R_only - unifSingAntData(t).FR_R_only(zeroI);
    analyzedData.pullL(t) = mean(FR_L(1:zeroI-1));
    analyzedData.pullR(t) = mean(FR_R(1:zeroI-1));
    analyzedData.pushL(t) = mean(FR_L(zeroI+1:end));
    analyzedData.pushR(t) = mean(FR_R(zeroI+1:end));
    
    analyzedData.pullX(t) = mean(unifSingAntData(t).aL_only(1:zeroI-1));
    analyzedData.pushX(t) = mean(unifSingAntData(t).aL_only(zeroI+1:end));
end

figure; hold on
a = gca;
a.ColorOrder = cat(1, a.ColorOrder, [38 34 97; 235 0 139; 0 104 56]./255);
plot([analyzedData.pullX; zeros(1,size(T,1)); analyzedData.pushX],...
     [analyzedData.pullL; zeros(1,size(T,1)); analyzedData.pushL], '-')
 ylim([-5 5])
export_fig(fullfile(scatfold, sprintf('allFlies_singleAnt_avgPULLavgPUSH_LEFTmoving.pdf')))

figure; hold on
a = gca;
a.ColorOrder = cat(1, a.ColorOrder, [38 34 97; 235 0 139; 0 104 56]./255);
plot([analyzedData.pullX; zeros(1,size(T,1)); analyzedData.pushX],...
     [analyzedData.pullR; zeros(1,size(T,1)); analyzedData.pushR], '-')
ylim([-5 5])
legend(num2str(T.flynum))
legend('Location', 'Best')
export_fig(fullfile(scatfold, sprintf('allFlies_singleAnt_avgPULLavgPUSH_RIGHTmoving.pdf')))



%% check you can replot scatter and interpolated 2d map with this unified data
clearvars -except unifSingAntData data T singleAntData scatfold
for t = 1:size(T,1)
    aR = data(t).aR;
    aL = data(t).aL;
    directions = data(t).dir;
    intensities = data(t).int;
    meanFR = data(t).FR;   
    
    if ~ismember(t, 8:10)
        zero_joint = 1;
        
        ixL = find(aR==0);
        ixL = ixL([2:4,1,5:7]);
        aL_ixL = aL(ixL);
        
        ixR = find(aL==0);
        ixR = ixR([2:4,1,5:7]);
        aR_ixR = aR(ixR);   % xaxis - intervals are not linear, plus values might not be sorted? - but they are
        
    else % 223 233 flies
        zero_joint = 25;
        
        ixL = find(aR==0);
        ixR = find(aL==0);
        aL_ixL = aL(aR==0); %these are already sorted!
        aR_ixR = aR(aL==0);  %these are already sorted!
        
    end
    
%% reg scatter    
%     figure; hold on; axis image
%     xlabel(' Left antenna (ipsi) - angular displ')
%     ylabel(' Right antenna (contra) - angular displ')
%     scatter(aL, aR, 200, meanFR-meanFR(zero_joint), 'filled', 'MarkerEdgeColor', 'k')
%     title(sprintf('Vm (N = %d)\n fly %d', T.NtrialsIncluded(t), T.flynum(t)), 'Interpreter', 'none')
%     colormap(bluewhitered(256)), cb = colorbar;
%     cb.TickLabels = num2str((cb.Ticks+meanFR(zero_joint))', '%2.1f');
%     cb.Label.String = '(Hz)';
%     xlim([-15 15])
%     ylim([-15 15])
%     set(gca, 'TickDir', 'out')    % OK!
%     
%% calculate scatter and interpolation:   

    % lin sum scatter
    [Lx,Ry] = meshgrid(aL_ixL, aR_ixR);
    [FRx,FRy] = meshgrid(meanFR(ixL), meanFR(ixR));
    Z = (FRx+FRy)/2 -meanFR(zero_joint);
    
    % lin sum interpolate
    F = scatteredInterpolant(Lx(:),Ry(:),Z(:));
    % F.Method = 'linear';
    % ExtrapolationMethod = 'nearest';
    F.Method = 'natural';
    ExtrapolationMethod = 'nearest';
    
    [Xq,Yq] = meshgrid(min(aL_ixL):0.25:max(aL_ixL));
    Vq = F(Xq, Yq);
    
    
    
    
    % original scatter and interpolate -- for difference map
%     aL;
%     aR;
%     meanFR-meanFR(zero_joint);
    V_data = meanFR'-meanFR(zero_joint);
    F_data = scatteredInterpolant(aL(:),aR(:),V_data);
    % F.Method = 'linear';
    % ExtrapolationMethod = 'nearest';
    F_data.Method = 'natural';
    ExtrapolationMethod = 'nearest';
    
    [Xqd,Yqd] = meshgrid(min(aL):0.25:max(aL));
    Vqd = F(Xqd, Yqd);
    
    
    
    
    XLIMs = [max([min(Xq(:)), min(Xqd(:))]), min([max(Xq(:)), max(Xqd(:))])];
    YLIMs = [max([min(Yq(:)), min(Yqd(:))]), min([max(Yq(:)), max(Yqd(:))])];
    
    % normalize:
    Nqd = Vqd + meanFR(zero_joint);
    Nq  =  Vq + meanFR(zero_joint);
    
    
    
%% plot scatter    
% %     figure; hold on; axis image
% %     mesh(Lx,Ry,Z)
% %     colormap(bluewhitered(256)), cb = colorbar;
% %     cb.TickLabels = num2str((cb.Ticks+meanFR(zero_joint))', '%2.1f');
% %     xlabel('L')
% %     ylabel('R')
% %     view(0, 90)
% %     

% already done:
    figure; hold on; axis image
    xlabel('L')
    ylabel('R')
    scatter(Lx(:), Ry(:), 200, Z(:), 'filled', 'MarkerEdgeColor', 'k')
    colormap(bluewhitered(256)), cb = colorbar;
    cb.TickLabels = num2str((cb.Ticks+meanFR(zero_joint))', '%2.1f');
    xlim([-15 15])
    ylim([-15 15])
    set(gca, 'TickDir', 'out') 
    title(sprintf('linear sum of FR of either antenna (divided by 2)\nfly %d', T.flynum(t)), 'Interpreter', 'none')
    export_fig(fullfile(scatfold, sprintf('linSum_scatter_%d.pdf', T.flynum(t))))
    
    
%% plot interpolation
    figure;
    
    subplot(1,2,1)
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
    cb.TickLabels = num2str((cb.Ticks+meanFR(zero_joint))', '%2.1f');
    cb.Label.String = '(Hz)';
    
    title(sprintf('interp from single antenna data - fly %d',T.flynum(t)), 'Interpreter', 'none')
% export_fig(fullfile(scatfold, sprintf('interp_linSum_scatter_%d.pdf', T.flynum(t))))
    
    subplot(1,2,2)
    V = meanFR'-meanFR(zero_joint);
    F = scatteredInterpolant(aL(:),aR(:),V);
    % F.Method = 'linear';
    % ExtrapolationMethod = 'nearest';
    F.Method = 'natural';
    ExtrapolationMethod = 'nearest';
    
    [Xq,Yq] = meshgrid(min(aL):0.25:max(aL));
    Vq = F(Xq, Yq);
    

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
    cb.TickLabels = num2str((cb.Ticks+meanFR(zero_joint))', '%2.1f');
    cb.Label.String = '(Hz)';
    
    title('interpolation of joint-antenna recorded data', 'Interpreter', 'none')
export_fig(fullfile(scatfold, sprintf('interp_2compare_%d.pdf', T.flynum(t))))
end











