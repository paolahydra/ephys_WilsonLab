%% table makng in 2017b
nCells = 9;  %abundanda.... or so I thought...

flynum = nan(nCells, 1);
cellfolder = cell(nCells, 1);
genotype = cell(nCells, 1);
hasFill = nan(nCells, 1);
metadatafiles = cell(nCells, 1);
SPIKEDET_STDMULTIPL = cell(nCells, 1);
NtrialsIncluded = nan(nCells, 1);
linkToDataFile = cell(nCells, 1);
notes = cell(nCells, 1);

meanFR_I30 = nan(nCells, 11);
meanFR_I60 = nan(nCells, 28);
meanFR_I120 = nan(nCells, 28);

meanVm_I30 = nan(nCells, 11);
meanVm_I60 = nan(nCells, 28);
meanVm_I120 = nan(nCells, 28);

directions_I30 = nan(nCells, 11);
directions_I60 = nan(nCells, 28);
directions_I120 = nan(nCells, 28);




T = table(  flynum, ...
            cellfolder,...
            genotype,...
            hasFill, ...
            metadatafiles,...
            SPIKEDET_STDMULTIPL, ...
            NtrialsIncluded, ...
            linkToDataFile, ...
            notes, ...
            meanFR_I30, meanFR_I60, meanFR_I120, ...
            meanVm_I30, meanVm_I60, meanVm_I120, ...
            directions_I30, directions_I60, directions_I120);
        
save('/Users/galileo/Dropbox (HMS)/p2/dataTable_p2.mat', 'T');

%%

load('/Users/galileo/Dropbox (HMS)/p2/dataTable_p2.mat');
% data folders
i = 1;
T.flynum(i)      = 292;
T.cellfolder{i}  = '/Users/galileo/Dropbox (HMS)/p2/fly292_PP/fly292_cell01'; 
T.genotype{i}    = 'R31A12';
T.hasFill(i)     = 0;
T.metadatafiles{i}          = {'all'};
T.SPIKEDET_STDMULTIPL{i}    = 2.5;
T.NtrialsIncluded(i)        = 40;


i = i+1;
T.flynum(i) = 290;
T.cellfolder{i} = '/Users/galileo/Dropbox (HMS)/p2/fly290_PP/fly290_cell01'; 
T.genotype{i}    = 'R31A12';
T.hasFill(i)     = 0;
T.metadatafiles{i} = {'metadata_20190201T173529.mat'};
T.SPIKEDET_STDMULTIPL{i} = 2.0;
T.NtrialsIncluded(i) = 11;



i = i+1;
T.flynum(i) = 267;
T.cellfolder{i} = '/Users/galileo/Dropbox (HMS)/p2/fly267_PP/fly267_cell01';
T.genotype{i}    = 'VT';
T.hasFill(i)     = 1;
T.metadatafiles{i} = {'all'};
T.SPIKEDET_STDMULTIPL{i} = 3.5;
T.NtrialsIncluded(i) = 14;


i = i+1;
T.flynum(i) = 264;
T.cellfolder{i} = '/Users/galileo/Dropbox (HMS)/p2/fly264_PP/fly264_cell01';
T.genotype{i}    = 'VT';
T.hasFill(i)     = 0; %did not fill
T.metadatafiles{i} = {'metadata_20181228T123901.mat'; 'metadata_20181228T130720.mat'};
T.SPIKEDET_STDMULTIPL{i} = [2.5, 2];
T.NtrialsIncluded(i) = 10;
T.notes{i} = {'cell did not fill, recording did not last very long'; ...
    'probes were repositioned and recalibrated between run 1 and run 2, because noted some instability of attachment of both antennae during run 1';
    'during run 2 the cell is much more depolarized than during run 1'};


i = i+1;
T.flynum(i) = 263;
T.cellfolder{i} = '/Users/galileo/Dropbox (HMS)/p2/fly263_PP/fly263_cell01';
T.genotype{i}    = 'VT';
T.hasFill(i)     = 1; 
T.metadatafiles{i} = {'all'}; %run 1: blocks 2:4; run 2: blocks 1:4
T.SPIKEDET_STDMULTIPL{i} = [2.5,2];
T.NtrialsIncluded(i) = 16;



i = i+1;
T.flynum(i) = 279;
T.cellfolder{i} = '/Users/galileo/Dropbox (HMS)/p2/fly279_PP/fly279_cell01';
T.genotype{i}    = 'VT';
T.hasFill(i)     = 0; %?
T.metadatafiles{i} = {'metadata_20190114T132434.mat';'metadata_20190114T133405.mat';'metadata_20190114T134314.mat'}; %run 1: blocks 2:4; run 2: blocks 1:4
T.SPIKEDET_STDMULTIPL{i} = [3, 4.5, 6];  % interspike interval increase from 5 to 8 ms.
T.NtrialsIncluded(i) = 4 + 4 + 2; %run 3: only keep first block, hoping that both reps are good (second block is too depolarized)



i = i+1;
T.flynum(i) = 283; 
T.cellfolder{i} = '/Users/galileo/Dropbox (HMS)/p2/fly298_PP/fly298_cell01';
T.genotype{i}    = 'VT';
T.hasFill(i)     = 1; % and it's anterior
T.metadatafiles{i} = {'all'}; 
T.SPIKEDET_STDMULTIPL{i} = [2.5, 2, 2];  % NOT positive about thresholds. interspike interval 8 ms.
T.NtrialsIncluded(i) = 23; 


i = 8;
T.flynum(i) = 223; 
T.cellfolder{i} = '/Users/galileo/Dropbox (HMS)/p2/fly223_PP/fly223_cell01';
T.genotype{i}    = 'VT';
T.hasFill(i)     = 0; %it was lost somewhere
T.metadatafiles{i} = {'metadata_20180509T193643.mat'; 'metadata_20180509T200130.mat'}; % of 193643 keeping all 7 blocks (of 3x), despite some of the very last stimuli were with piezo detatching from antennae
T.SPIKEDET_STDMULTIPL{i} = [3.5];  %  interspike interval 8 ms.
T.NtrialsIncluded(i) = 21+6; 
T.notes{i} = {'for 193643 estimated: armL=320; armR=350;  '; ...
    '200130 had moved posterior (R) piezo more proximal -> armR = 320'};

i = 9;
T.flynum(i) = 2331; %fly 233 cell 1
T.cellfolder{i} = '/Users/galileo/Dropbox (HMS)/p2/fly233_PP/fly233_cell01';
T.genotype{i}    = 'VT';
T.hasFill(i)     = 0; %it was lost somewhere
T.metadatafiles{i} = {'all'}; % keeping 4 blocks out of 7. Cell is very depolaryzed and dying at the end 
T.SPIKEDET_STDMULTIPL{i} = [1.25];  %  interspike interval 10 ms.
T.NtrialsIncluded(i) = 12; 
T.notes{9} = {'arms are shorter than 223, but vertical alignment much worse, esp for anterior/L antenna'};


i = 10;
T.flynum(i) = 2332; %fly 233 cell 2
T.cellfolder{i} = '/Users/galileo/Dropbox (HMS)/p2/fly233_PP/fly233_cell02';
T.genotype{i}    = 'VT';
T.hasFill(i)     = 0; %it was lost somewhere
T.metadatafiles{i} = {'all'}; % of 193643 keeping all 7 blocks (of 3x), despite some of the very last stimuli were with piezo detatching from antennae
T.SPIKEDET_STDMULTIPL{i} = [0.9, 0.9];  %  interspike interval 10 ms. - both few false positive, and more false negative in bursts
T.NtrialsIncluded(i) = 6+21; 
% T.notes{i} = {};


% i = i+1;
% T.flynum(i) = 298; %HUGE FLUCTUATIONS IN VM (ESPECIALLY RUN 1)
% T.cellfolder{i} = '/Users/galileo/Dropbox (HMS)/p2/fly298_PP/fly298_cell01';
% T.genotype{i}    = 'VT';
% T.hasFill(i)     = 0; %?
% T.metadatafiles{i} = {'metadata_20190211T114642.mat'; 'metadata_20190211T115222.mat'}; 
% T.SPIKEDET_STDMULTIPL{i} = [4];  % interspike interval increase from 5 to 8 ms.
% T.NtrialsIncluded(i) = 2 + 8; 


save('/Users/galileo/Dropbox (HMS)/p2/dataTable_p2.mat', 'T', '-append');
