figure; hold on
for t = [1,5]
%     plot(T.directions_I120(t,:), T.meanFR_I120(t,:),'-o', 'MarkerFaceColor', 'auto', 'MarkerSize', 3)
    plot(T.directions_I60(t,:), T.meanFR_I60(t,:),'-o', 'MarkerFaceColor', 'auto', 'MarkerSize', 3)
%     plot(T.directions_I30(t,:), T.meanFR_I30(t,:),'-o', 'MarkerFaceColor', 'auto', 'MarkerSize', 3)
end
legend(num2str(T.flynum(1:6)))
a = gca;
a.TickDir = 'out';
ylabel('baseline subtracted firing rate')
xlabel('wind direction (deg)')
export_fig('dtAVG_directionTuning_I120_2cells.pdf')

figure; hold on
for t = 1:6
    plot(T.directions_I60(t,:), T.meanFR_I60(t,:),'-o', 'MarkerFaceColor', 'auto', 'MarkerSize', 3)
end

figure; hold on
for t = 1:6
    plot(T.directions_I30(t,:), T.meanFR_I30(t,:),'-o', 'MarkerFaceColor', 'auto', 'MarkerSize', 3)
end




% export_fig('/Users/galileo/Dropbox (HMS)/OkuboPatellaWilson ms/pp_bits/AVG_dirTuning_firingRateRestSubtr_fly292_cell01.pdf')