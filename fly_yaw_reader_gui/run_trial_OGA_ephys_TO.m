function [ trial_data, trial_time ] = run_trial_OGA_ephys_TO( trial_idx, task, run_obj, scanimage_client, trial_core_name )
%%% Tatsuo Okubo
%%% fictive odor-gated anemotaxis with ephys
%%% 2017/01/04

disp(['About to start trial task: ' task]);

Ch = rig_params_TO;

% Setup data structures for read / write on the daq board
s = daq.createSession('ni');

% This channel is for external triggering of scanimage 5.1
s.addDigitalChannel(Ch.DevID, 'port0/line0', 'OutputOnly');

% These are for inputs: motion sensor 1 x,y; motion sensor 2 x,y; frame
% clock; stim left; stim right;
ai_channels_used = [0:9];
aI = s.addAnalogInputChannel(Ch.DevID, ai_channels_used, 'Voltage');
for i=1:length(ai_channels_used)
    aI(i).InputType = 'SingleEnded';
end

% This is the stim control: stim left, stim right
s.addAnalogOutputChannel(Ch.DevID, 0:1, 'Voltage'); % current injection, shutter for wide field illumination
s.addDigitalChannel(Ch.DevID, ['port0/line1:4'], 'OutputOnly'); % additional channel for the valves (air, L/C/R)

settings = sensor_settings;

SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
total_duration = run_obj.pre_stim_t + run_obj.stim_t + run_obj.post_stim_t;

% initialize the output vectors to zero
zero_stim = zeros(SAMPLING_RATE*total_duration,1);
stim_LED = zeros(SAMPLING_RATE*total_duration,1);
stim_wind = zeros(SAMPLING_RATE*total_duration,1);
stim_pinch = zeros(SAMPLING_RATE*total_duration,1);

% set the stimulation pulse
begin_idx = run_obj.pre_stim_t * SAMPLING_RATE;
end_idx = (run_obj.pre_stim_t+run_obj.stim_t) * SAMPLING_RATE;
Margin.pre = 2; % offset between trial onset and wind onset [s]
Margin.post = 1; % offset between wind offset and trial offset [s]
begin_wind = Margin.pre * SAMPLING_RATE;
end_wind = (total_duration - Margin.post) * SAMPLING_RATE;
stim_LED(begin_idx:end_idx) = 5.0; % analog output for LED
stim_wind(begin_wind:end_wind) = 1; % for master valve that controls the air flow
stim_pinch(2:(end-1))=1; % pinch valves: high or low throughout the trial

output_data = [];

if( strcmp(task, 'OdorLeftWind') == 1 )    
    output_data = [stim_LED stim_LED zero_stim stim_wind zero_stim zero_stim stim_pinch];
elseif( strcmp(task, 'OdorRightWind') == 1 )
    output_data = [stim_LED stim_LED zero_stim stim_wind stim_pinch zero_stim zero_stim];
elseif( strcmp(task, 'OdorNoWind') == 1 )
    output_data = [stim_LED stim_LED zero_stim zero_stim zero_stim zero_stim zero_stim];
else
    disp(['ERROR: Task: ' task ' is not recognized.']);
end   

queueOutputData(s, output_data);

% Delay starting the aquisition for a second to ensure that scanimage is
% ready
pause(1.0);

[trial_data, trial_time] = s.startForeground();

release(s);
end