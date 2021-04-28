function varargout = setDCoffset_2Piezos()
%set DAC output channels and check piezo's DC offset

HWsettings = sensor_settings_PP; % to set AO channels consistently

daqreset;
s = daq.createSession('ni');
s.Rate = 1e4;

ai = s.addAnalogInputChannel(HWsettings.DAQdevName, ...
            [HWsettings.AIch.PIEZO_SENSOR_LEFT, HWsettings.AIch.PIEZO_SENSOR_RIGHT] , 'Voltage'); 
for ch = 1:2
    ai(ch).Range = [-10 10];
    ai(ch).TerminalConfig = 'SingleEnded';
end

%% tune DC offset
DC_offset = [0,0];
disp('turn SERVO ON')
while sum(DC_offset < [4.9, 4.9] | DC_offset > [5.1, 5.1]) 
    disp('tune piezos DCoffset between 4.9 and 5.1 and click any key:')
    pause
    read_in = startForeground(s);
    DC_offset = mean(read_in);
    fprintf(' LEFT  (anterior) DC:\t%1.2f\n', DC_offset(1))
    fprintf('RIGHT (posterior) DC:\t%1.2f\n', DC_offset(2))
end

%%
if nargout
    varargout{1} = DC_offset;
end
release(s)

end

