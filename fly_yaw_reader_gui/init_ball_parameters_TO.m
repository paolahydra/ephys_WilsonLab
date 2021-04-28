function [  ] = init_ball_parameters_TO(experiment_dir)
% init_ball_parameters: 
% WARNING: THE BALL MUST NOT BE MOVING
settings = sensor_settings;

session_obj = daq.createSession('ni');
Ch = rig_params_TO;

% These are for inputs: motion sensor 1 x,y; motion sensor 2 x,y; frame
% clock; stim left; stim right;
ai_channels_used = [0:15];
aI = session_obj.addAnalogInputChannel( Ch.DevID, ai_channels_used, 'Voltage' );
for i=1:length(ai_channels_used)
    aI(i).InputType = 'SingleEnded';
end

[trial_data, trial_time] = session_obj.startForeground();

zero_mean_per_channel = [];
zero_std_per_channel = [];

[ zero_mean_per_channel(1), zero_std_per_channel(1) ] = get_zero_velocity_for_channel( trial_data( :, settings.sensor_1_DX_DAQ_AI ) );
[ zero_mean_per_channel(2), zero_std_per_channel(2) ] = get_zero_velocity_for_channel( trial_data( :, settings.sensor_1_DY_DAQ_AI ) );
[ zero_mean_per_channel(3), zero_std_per_channel(3) ] = get_zero_velocity_for_channel( trial_data( :, settings.sensor_2_DX_DAQ_AI ) );
[ zero_mean_per_channel(4), zero_std_per_channel(4) ] = get_zero_velocity_for_channel( trial_data( :, settings.sensor_2_DY_DAQ_AI ) );


filename = settings.zero_params_filename;
save([experiment_dir '\' filename], 'zero_mean_per_channel', 'zero_std_per_channel');

release( session_obj );

end

