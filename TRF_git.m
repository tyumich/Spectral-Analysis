%% Loads fieldtrip toolbox

addpath('/MATLAB/fieldtrip-20170505');
ft_defaults
%% Setup subject number & Trigger

dataset='/MandarinRC_01.eeg';

if ~exist('proc')|| ~isfield(proc, 'subject')
    [tokens, pos] = regexp(dataset, 'R([\d]{4})', 'tokens');
    sidx = pos(1);

    proc.subject = dataset(sidx:sidx+4);
    proc.dataset = dataset;
end 

%trigger numbers

Conds = {'S 13' 'S 14' 'S 15' 'S 16' 'S 17' 'S 18' 'S 19' 'S 20' 'S 21' 'S 22' 'S 23'
         'S 33' 'S 34' 'S 35' 'S 36' 'S 37' 'S 38' 'S 39' 'S 40' 'S 41' 'S 42' 'S 43'
         'S 51' 'S 52' 'S 53' 'S 54' 'S 55' 'S 56' 'S 57' 'S 58' 'S 59' 'S 60' 'S 61'
         'S 71' 'S 72' 'S 73' 'S 74' 'S 75' 'S 76' 'S 77' 'S 78' 'S 79' 'S 80' 'S 81'};  
     
% Epoch & Load
% do each condition separately to manually add trialinfo

raw = {};

cfg = [];
cfg.dataset = dataset;
cfg.trialdef.eventtype = 'Stimulus';
cfg.trialdef.eventvalue = Conds;
cfg.trialdef.prestim = 0.3;   
cfg.trialdef.poststim = 1.0;    
cfg = ft_definetrial(cfg);

ntrials = length(cfg.trl);

% Load & filter
cfg.dataset = dataset;
cfg.channel = {'all', '-VEOG', '-AUD', '-OPTO'};
cfg.padding = 1; % 1 sec of padding for filters
cfg.implicitref = '29';
cfg.reref = 'yes';
cfg.refchannel = {'25', '29'};
cfg.hpfilter = 'yes';
cfg.hpfreq = 0.1;
cfg.hpfiltord = 3; 
cfg.dftfilter = 'yes';
cfg.dftfreq = [60 120 180];
%cfg.precision = 'single';

raw = ft_preprocessing(cfg);
%% View Data

cfg = [];
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, raw);
%raw2 = ft_databrowser(cfg, raw);

cfg        = [];
cfg.method = 'channel';
ft_rejectvisual(cfg, raw)
%% Make trialinfo and combine data

data_all = raw;
proc.trl = raw.cfg.trl;
%% Mark/remove high impedence chans
                                                                                                                                                                                                                                   
if ~exist('proc') || ~isfield(proc, 'impedence')
    [proc.impedence.bads proc.impedence.imps proc.impedence.labels] = get_high_impedence(dataset, 25);
    
    picks = setdiff(data_all.label, proc.impedence.bads);
    cfg = [];
    cfg.channel = picks;
    data_all = ft_selectdata(cfg, data_all);
    
end
%% Manually reject artifacts - initial sweep
% data_all -> data_rej1
% only remove really exceptionally bad chans (e.g. zeros)

if ~exist('proc') || ~isfield(proc, 'first_artfctdef')
    dummy                = ft_rejectvisual([], data_all);
    proc.first_artfctdef = dummy.cfg.artfctdef;
    proc.first_picks = dummy.label;
    clear dummy
else
    proc.first_picks = data_all.label;
end

cfg = [];
cfg.artfctdef = proc.first_artfctdef;
data_rej1 = ft_rejectartifact(cfg, data_all);

cfg = [];
cfg.channel = proc.first_picks;
data_rej1 = ft_selectdata(cfg, data_rej1);
%% ICA
% data_rej1 -> data_ica
% unmixing matrix computed over a downsampled dataset

artifact = proc.first_artfctdef.summary.artifact;
picks = proc.first_picks;

if ~exist('proc') || ~isfield(proc, 'ica')
    [proc.ica.unmixing, proc.ica.topolabel, proc.ica.rej_comp, proc.ica.comments] = ...
        get_alice2_eeg_ica(dataset, artifact, picks);
        %get_mandarin_ica(data_rej1, proc.first_artfctdef, proc.first_picks);
end

% unmix the lightly cleaned data...
cfg = [];
cfg.unmixing = proc.ica.unmixing;
cfg.topolabel = proc.ica.topolabel;
comp = ft_componentanalysis(cfg, data_rej1);

% ...then reject components
cfg = [];
cfg.component = proc.ica.rej_comp; % Excluded Components
data_ica = ft_rejectcomponent(cfg, comp, data_rej1);

clear comp
%% Manual trial rejections - final sweep
% data_ica -> data_rej2
% mark bad chans here

if ~exist('proc') || ~isfield(proc, 'second_artfctdef')
    dummy                                   = ft_rejectvisual([], data_ica);

    proc.second_artfctdef = dummy.cfg.artfctdef;
    proc.second_picks = dummy.label;
    
    clear dummy
end

cfg = [];
cfg.artfctdef = proc.second_artfctdef;
data_rej2 = ft_rejectartifact(cfg, data_ica);

cfg = [];
cfg.channel = proc.second_picks;
data_rej2 = ft_selectdata(cfg, data_rej2);

proc.badchans = [setdiff(data_all.label, proc.second_picks); proc.impedence.bads];
rejected_chans = setdiff(data_all.label, proc.second_picks);
proc.badchans = [rejected_chans(:); proc.impedence.bads(:)];

%   tracks all rejected chans + high impedences
proc.numtrialrej = length(data_all.trial) - length(data_rej2.trial);
%% Bad channels are replaced with nearest neighbour interpolation
% data_rej2 stays data_rej2

% track rank for further data decomposition (e.g. ica, beamforming &c.)
proc.rank = length(data_all.label) - ... % bad impedences already removed
            length(proc.badchans) - ... % rejected chans + bad impedences
            length(proc.ica.rej_comp);
%            length(proc.impedence.bads); % don't double-count bad impedences!

if ~isempty(proc.badchans);
    cfg         = [];
    cfg.method  = 'template';
    cfg.template                        = 'easycapM10-acti61_neighb.mat';
    neighbs = ft_prepare_neighbours(cfg);
    
    cfg         = [];
    cfg.method  = 'spline';
    cfg.badchannel = proc.badchans;
    cfg.neighbours = neighbs;
    cfg.elecfile                        = 'easycapM10-acti61_elec.sfp';
    data_rej2 = ft_channelrepair(cfg, data_rej2);
end
%% Separate into conditions
% data_rej2 -> dat

dat = {};
condition_integers = unique(data_rej2.trialinfo);

for c = 1:length(condition_integers)   
    cfg = [];
    cfg.trials = find(data_rej2.trialinfo == condition_integers(c));
    dat{c} = ft_selectdata(cfg, data_rej2);
end
%% TFRs 

% 5-80 Hz at 2 Hz intervals; 20 ms windows

tfr = {};
for c = 1:length(dat)
    cfg                 = [];
    cfg.method          = 'wavelet';
    cfg.output          = 'pow';
    cfg.foi             = 1:2:80;
    cfg.toi             = -2:.02:2;
    cfg.width           = 7;
    
    cfg.pad = 'nextpow2';
    cfg.padtype = 'zero';
    
    tfr{c} = ft_freqanalysis(cfg, dat{c});

    cfg = [];
    cfg.baseline = [-.3 -.1];
    cfg.baselinetype = 'relchange';
    cfg.parameter = 'powspctrm';

    tfr{c} = ft_freqbaseline(cfg, tfr{c});
end

% Visualization/Plot the result
cfg = [];
cfg.showlabels   = 'yes';
cfg.layout       = 'easycapM10-acti61.lay';
figure; ft_multiplotTFR(cfg, tfr{2}); %set trigger

cfg = [];
cfg.marker       = 'on';
cfg.layout       = 'easycapM10-acti61.lay';
cfg.channel = 'all';
% cfg.channel = {'32', '59', '60', '61', '20'}; %Left Anterior
cfg.interactive  = 'no';
cfg.ylim         = [5 30]; %set freq
cfg.zlim         = [0 20]; %set powspctrm
figure;
subplot(211);ft_singleplotTFR(cfg, tfr{7}); title('MV: SS'); %set trigger
subplot(212);ft_singleplotTFR(cfg, tfr{18}); title('MV: SO'); %set trigger

cfg = [];
% cfg.xlim         = [.1 .3]; %set time
% cfg.zlim         = [0 10]; %set powspctrm
% cfg.ylim         = [3 8]; %set freq
cfg.marker       = 'on';
cfg.layout       = 'easycapM10-acti61.lay';
figure;
subplot(211);ft_topoplotTFR(cfg, tfr{7}); title('MV: SS'); %set trigger
subplot(212);ft_topoplotTFR(cfg, tfr{18}); title('MV: SO'); %set trigger

%Part II: Spectral analysis on EEG resting state data
%ft_redefinetrial to cut data into non-overlapping segments of various lengths (1 sec, 2 secs and 4 secs)
%compute power spectrum of all data segments and average them

cfg1 = [];
cfg1.length  = 1; %set segment length
cfg1.overlap = 0;
base_rpt1    = ft_redefinetrial(cfg1, dattl{2}); %set trigger

cfg2 = [];
cfg2.output  = 'pow';
cfg2.channel = 'all';
cfg2.method  = 'mtmfft';
cfg2.taper   = 'boxcar';
cfg2.foi     = 0.5:1:45; % 1/cfg1.length  = 1;
base_freq1   = ft_freqanalysis(cfg2, base_rpt1);

figure;
hold on;
plot(base_freq1.freq, base_freq1.powspctrm(50,:)) %set channel
xlabel('Frequency (Hz)');
ylabel('absolute power (uV^2)');
% plot(base_freq1.freq, base_freq1.powspctrm(61,:))
% plot(base_freq2.freq, base_freq2.powspctrm(61,:))
% plot(base_freq4.freq, base_freq4.powspctrm(61,:))
% legend('1 sec window','2 sec window','4 sec window')

figure;
hold on;
plot(tfr{2}.freq, tfr{2}.powspctrm(61,:)) %set channel
xlabel('Frequency (Hz)');
ylabel('absolute power (uV^2)');

%% Save Data

save('R0000_01.mat','dat','proc')
