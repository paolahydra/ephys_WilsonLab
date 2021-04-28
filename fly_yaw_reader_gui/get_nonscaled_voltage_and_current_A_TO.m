function [ current, voltage ] = get_nonscaled_voltage_and_current_A_TO( trial_data )
%%% get current and voltage directly (i.e. not through scaled output)
%%% Tatsuo Okubo
%%% 2017/04/13

settings = sensor_settings;

mode_meta = mean(trial_data(:,settings.MODE_A_DAQ_AI));
mode = get_mode_for_meta_TO( mode_meta );

current_nonscaled = trial_data(:,settings.CURRENT_NON_SCALED_A_DAQ_AI); % [V] 
voltage_nonscaled = trial_data(:,settings.VOLTAGE_NON_SCALED_A_DAQ_AI)/10; % [V], note 10Vm
gain_meta = mean(trial_data(:,settings.GAIN_A_DAQ_AI));

current = current_nonscaled*1000; % [pA], based on the rear panel settings
voltage = voltage_nonscaled*1000; % [mV]