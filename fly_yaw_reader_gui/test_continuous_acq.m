function [s] = view_fly_trajectory()

s = daq.createSession('ni');

ai_channels_used = [0:3];

aI = s.addAnalogInputChannel('Dev1', ai_channels_used, 'Voltage');

for i=1:length(ai_channels_used)
    aI(i).InputType = 'SingleEnded';
end

settings = sensor_settings;
s.Rate = settings.sampRate;
s.NotifyWhenDataAvailableExceeds = 2000;
myfig = figure();
lh = addlistener(s,'DataAvailable', @(src,event)display_fly_trajectory(src,event,myfig));
s.IsContinuous = true;

global last_pos_x;
global last_pos_y;
global last_theta;
global iter_count;
    
last_pos_x = 0;
last_pos_y = 0;
last_theta = 0;
iter_count = 0;

s.startBackground()

end