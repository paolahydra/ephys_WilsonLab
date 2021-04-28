function HWsettings = sensor_settings_PP_olderFlies
% this collects all hardwired settings that are specific to a given rig
% structure changed on 4/19/2018 - branching away from TO's structure
%
% simply returns the AI channel number relative to specific inputs.
% channel numbers are in 0+ coordinates
% matlab data is in 1+ coordinates.
% this will be used to retrieve matlaba data, therefore adding 1 to the
% channel number.
 
 
 
HWsettings.DAQdevName  = 'Dev1'; %just the name
 
 
 
%% AI: put your wiring number here (channel numbers)
HWsettings.AIch.MODE_A_DAQ_AI                 = 0;
HWsettings.AIch.VOLTAGE_NON_SCALED_A_DAQ_AI   = 1;
HWsettings.AIch.SCALED_OUT_A_DAQ_AI           = 2;
HWsettings.AIch.GAIN_A_DAQ_AI                 = 3;
HWsettings.AIch.FREQUENCY_A_DAQ_AI            = 4;        %check that this name consistently used, if ever used
HWsettings.AIch.CURRENT_NON_SCALED_A_DAQ_AI   = 5;
 
HWsettings.AIch.PIEZO_SENSOR_LEFT             = 6;
HWsettings.AIch.PIEZO_SENSOR_RIGHT            = 7;
 
% HWsettings.AIch.SHUTTER                       = 8;        %what is this for? 
 
 
 
%% AO:
HWsettings.AOch.PIEZO_OUTPUT_LEFT  = 0; % analogous output channel
HWsettings.AOch.PIEZO_OUTPUT_RIGHT = 1; % analogous output channel
 
 
 
%% data retrieving settings (matlab 1- based indexing)
HWsettings.DATAch.MODE_A_DAQ_AI               = HWsettings.AIch.MODE_A_DAQ_AI               +1;
HWsettings.DATAch.GAIN_A_DAQ_AI               = HWsettings.AIch.GAIN_A_DAQ_AI               +1;
HWsettings.DATAch.FREQUENCY_A_DAQ_AI          = HWsettings.AIch.FREQUENCY_A_DAQ_AI          +1;
HWsettings.DATAch.CURRENT_NON_SCALED_A_DAQ_AI = HWsettings.AIch.CURRENT_NON_SCALED_A_DAQ_AI +1;
HWsettings.DATAch.VOLTAGE_NON_SCALED_A_DAQ_AI = HWsettings.AIch.VOLTAGE_NON_SCALED_A_DAQ_AI +1;
HWsettings.DATAch.SCALED_OUT_A_DAQ_AI         = HWsettings.AIch.SCALED_OUT_A_DAQ_AI         +1;
 
HWsettings.DATAch.PIEZO_SENSOR_LEFT           = HWsettings.AIch.PIEZO_SENSOR_LEFT           +1;
HWsettings.DATAch.PIEZO_SENSOR_RIGHT          = HWsettings.AIch.PIEZO_SENSOR_RIGHT          +1;
% HWsettings.DATAch.SHUTTER                     = HWsettings.AIch.SHUTTER                     +1;
 
 
 
HWsettings.DATAsett.CURRENT_GAIN       = 10;
HWsettings.DATAsett.CURRENT_FREQUENCY  = 5;    % KHz
HWsettings.DATAsett.VOLTAGE_GAIN       = 10;
HWsettings.DATAsett.VOLTAGE_FREQUENCY  = 5;    % KHz
 
 
end