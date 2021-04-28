%% Test daq boards
settings = sensor_settings;
s = daq.createSession('ni');

aI = s.addAnalogInputChannel('Dev1', 0:3, 'Voltage');
s.Rate = settings.sampRate;
s.DurationInSeconds = 20;

for i=1:4
aI(i).InputType = 'SingleEnded';
end

[data, time] = s.startForeground();

figure;
for i=1:4
    subplot(4,1,i)
    plot(time, data(:,i));
    xlim([0 20]);
end

s.release()

%% Velocity when the ball is not moving

%2*std(data)
%mean(data)

two_std_per_channel = [0.0137, 0.0151, 0.0142, 0.0156];
zero_mean_voltage_per_channel = [2.4606, 2.4466, 2.4703, 2.4665];



%%
figure;

for i=1:4

subplot(4,1,i);    
raw_data = data(:,i);

rate = 2*(settings.cutoffFreq/settings.sampRate);
[kb, ka] = butter(2,rate);
smoothedData = filtfilt(kb, ka, raw_data);
dt = settings.sampRate/settings.sensorPollFreq;

smoothedData_downsampled = squeeze(mean(reshape(smoothedData, [dt, length(smoothedData)/dt])));
time_downsampled = squeeze(mean(reshape(time, [dt, length(time)/dt])));

smoothedData_downsampled_zeroed = smoothedData_downsampled - repmat(zero_mean_voltage_per_channel(i), [1 size(smoothedData_downsampled,2)]);

%smoothedData_downsampled_zeroed(find( (smoothedData_downsampled_zeroed < two_std_per_channel(i)) & (smoothedData_downsampled_zeroed > -1*two_std_per_channel(i)))) = 0.0;

hold on;
%plot(time, raw_data, 'color', 'b')
%plot(time, smoothedData, 'color', 'r')
%plot(time_downsampled, smoothedData_downsampled, 'color', 'g');
plot(time_downsampled, smoothedData_downsampled_zeroed, 'color', 'r');
ylim([-0.3, 0.3]);

end




if 0
voltsPerStep = (maxVal - minVal)/(settings.numInts - 1);
seq = round((smoothedData - minVal)./voltsPerStep);
%seq = (smoothedData - minVal)./voltsPerStep;
maxInt = settings.numInts -1;
seq(seq>maxInt) = maxInt;
seq(seq<0) = 0;
zeroVal = -1 + (settings.numInts + 1)/2;
seq = seq - zeroVal;

vel_in_mm = seq.*settings.mmPerCount.*settings.sensorPollFreq;
end
