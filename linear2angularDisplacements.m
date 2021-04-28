% calculate angular displacements
l = 200 ;        %um approximate length of arista arm at point of piezo attachment.

% displacements       = [-22,  -11,  -5,   0,   5,   11,   22]; % um
displacements       = [-28, -14,  -7,   0,   7,   14,  28]; % um




alphas = asind(displacements./l)
