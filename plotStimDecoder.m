figure;
hold on
% scatter(aL,aR)
axis square

text(aL, aR, num2str((1:84)'))
xlim([-15, 15])
ylim([-15, 15])

export_fig('stimIdx84_decoder.pdf')