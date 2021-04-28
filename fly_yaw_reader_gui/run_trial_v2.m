function [ trial_data, trial_time ] = run_trial_v2( trial_idx, task, run_obj, trial_core_name )

disp(['About to start trial task: ' task]);

% Setup data structures for read / write on the daq board
s = daq.createSession('ni');

s.addAnalogOutputChannel('Dev2', 0:1, 'Voltage');
%s.addDigitalChannel('Dev2', 'port0/line6', 'OutputOnly');
%s.addDigitalChannel('Dev2', 'port0/line7', 'OutputOnly');

% These are for inputs: motion sensor 1 x,y; motion sensor 2 x,y; frame
% clock; stim left; stim right;
ai_channels_used = [0:13];
aI = s.addAnalogInputChannel('Dev2', ai_channels_used, 'Voltage');
for i=1:length(ai_channels_used)
    aI(i).InputType = 'SingleEnded';
end

% This is the stim control: stim left, stim right
%

settings = sensor_settings;

SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
total_duration = run_obj.pre_stim_t + run_obj.stim_t + run_obj.post_stim_t;

zero_stim = zeros(SAMPLING_RATE*total_duration,1);
stim = zeros(SAMPLING_RATE*total_duration,1);

begin_idx = run_obj.pre_stim_t * SAMPLING_RATE;
end_idx = (run_obj.pre_stim_t+run_obj.stim_t) * SAMPLING_RATE;
stim(begin_idx:end_idx) = 5.0;

output_data = [];
if( strcmp(task, 'LeftOdor') == 1 )
    output_data = [ stim zero_stim ];
elseif( strcmp(task, 'RightOdor') == 1 )
    output_data = [ zero_stim stim ];
elseif( strcmp(task, 'BothOdor') == 1 )
    output_data = [ stim stim ];
elseif( strcmp(task, 'ExternalCommand') == 1 )
    
    CURRENT_INJECTION_IN_pA_PER_VOLT = 127.0; % This number assumes a 0.0637 voltage drop across the downstep resistors
    
    injection_current = run_obj.injection_current;
    
    external_command = zeros( SAMPLING_RATE*total_duration , 1 ); 
    
    external_command(begin_idx:end_idx) = injection_current ./ CURRENT_INJECTION_IN_pA_PER_VOLT;
    
    output_data = [ external_command zero_stim ];
    
elseif( strcmp(task, 'NaturalOdor') == 1 )
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is where olfactometer stim parameters are defined.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    PINCH_VALVE_OPEN_TIME = 1.0;
    begin_idx = PINCH_VALVE_OPEN_TIME * SAMPLING_RATE;
    
    pinch_valve_waveform = zeros(SAMPLING_RATE*total_duration,1);
    pinch_valve_waveform( begin_idx:end-1 ) = 5.0; % volts
    
    output_data = [pinch_valve_waveform stim ];
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
else
    disp(['ERROR: Task: ' task ' is not recognized.']);
end

queueOutputData(s, output_data);

[trial_data, trial_time] = s.startForeground();

release(s);
end

