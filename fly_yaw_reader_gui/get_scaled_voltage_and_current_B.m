function [ current, voltage ] = get_scaled_voltage_and_current_B( trial_data )

settings = sensor_settings;

mode_meta = mean(trial_data(:,settings.MODE_B_DAQ_AI));
mode = get_mode_for_meta( mode_meta );

current_nonscaled = (trial_data(:,settings.CURRENT_NON_SCALED_B_DAQ_AI)/10)*1000; % pA, gain of 10 from the back
voltage_nonscaled = (trial_data(:,settings.VOLTAGE_NON_SCALED_B_DAQ_AI)/10)*1000; % mV, gain of 10 from the back
scaledOut = trial_data(:,settings.SCALED_OUT_B_DAQ_AI);
gain_meta = mean(trial_data(:,settings.GAIN_B_DAQ_AI));

gain = get_gain_for_meta(gain_meta);

% Convention is mode == the value that ScaledOut represents
if( mode == 'V' )
    current = current_nonscaled;
    voltage = (scaledOut/gain) * 1000; % mV    
elseif( mode == 'I' )
    voltage = voltage_nonscaled;
    current = (scaledOut/gain) * 1000; % pA    
end

end

