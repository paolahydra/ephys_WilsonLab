function [ current, voltage, scaled, mode ] = get_scaled_voltage_and_current_PP( data )
%%% output is in [pA] for current, [mV] for voltage
%%% Tatsuo Okubo
%%% 2017/07/06

HWsettings = sensor_settings_PP;

mode_meta = mean(data(:,HWsettings.DATAch.MODE_A_DAQ_AI)); 
mode = get_mode_for_meta_TO( mode_meta ); % mode == the value that ScaledOut represents:
% mode = 'V'; % current clamp (note_PP: I0, normal or fast could be further disambiguated), scaledOutput is Vm
% mode = 'I'; % voltage clamp, scaledOutput is I <-> voltage clamp


gain_meta = mean(data(:,HWsettings.DATAch.GAIN_A_DAQ_AI));
[I_gain, V_gain] = get_gain_TO(gain_meta); %not sure why this read - because acq was not singleEndedd- careful with that
%


%%
current_nonscaled = data(:, HWsettings.DATAch.CURRENT_NON_SCALED_A_DAQ_AI); % raw recorded current [V] 
voltage_nonscaled = data(:, HWsettings.DATAch.VOLTAGE_NON_SCALED_A_DAQ_AI); % raw recorded voltage [V]
scaledOut = data(:, HWsettings.DATAch.SCALED_OUT_A_DAQ_AI); % scaled output [V]

%resistance is correctly calculated as long as the external filter gains for I and V are set
%the same. -PP However, in order to be faithful to the actual values, I am
%considering the scaling values in the external filter/amplifier as well
%(recorded manually -- try not to change them).

current = (current_nonscaled*1000) / HWsettings.DATAsett.CURRENT_GAIN; % [pA] based on 1 mV/pA conversion 
voltage = (voltage_nonscaled*1000/10) / HWsettings.DATAsett.VOLTAGE_GAIN ; % [mV], note that the output is 10Vm (PP- correct)


% Convention is mode == the value that ScaledOut represents
if strcmp(mode,'V')
    scaled = (scaledOut/V_gain) * 1000; % mV    
elseif strcmp(mode,'I')
    scaled = (scaledOut/I_gain) * 1000; % pA    
end


end