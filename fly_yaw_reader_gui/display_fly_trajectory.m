function  display_fly_trajectory( src, event, myfig )

figure(myfig);

settings = sensor_settings;

sensor1_x = event.Data(:,1);
sensor1_y = event.Data(:,2);
sensor2_x = event.Data(:,3);
sensor2_y = event.Data(:,4);

[ t, vel1_x ] = get_velocity_from_raw_input( sensor1_x, event.TimeStamps, settings.zero_mean_voltage_per_channel(1), settings.zero_mean_two_std_per_channel(1) );
[ t, vel1_y ] = get_velocity_from_raw_input( sensor1_y, event.TimeStamps, settings.zero_mean_voltage_per_channel(2), settings.zero_mean_two_std_per_channel(2) ); 
[ t, vel2_x ] = get_velocity_from_raw_input( sensor2_x, event.TimeStamps, settings.zero_mean_voltage_per_channel(3), settings.zero_mean_two_std_per_channel(3) );
[ t, vel2_y ] = get_velocity_from_raw_input( sensor2_y, event.TimeStamps, settings.zero_mean_voltage_per_channel(4), settings.zero_mean_two_std_per_channel(4) );

%vel_forward = -1*((vel1_y + vel2_y)*cos(deg2rad(45)));
%vel_side    = -1*((vel1_y - vel2_y)*sin(deg2rad(45)));
%vel_yaw     = -1*((vel1_x + vel2_x) ./ 2.0);
vel_forward = ((vel1_y + vel2_y)*cos(deg2rad(45)));
vel_side    = ((vel2_y - vel1_y)*sin(deg2rad(45)));
vel_yaw     = ((vel1_x + vel2_x) ./ 2.0);

global last_pos_x;
global last_pos_y;

global last_pos_x1;
global last_pos_y1;
global last_theta1;

global iter_count;

if(mod(iter_count, 10000) == 0)
    last_pos_x = 0;
    last_pos_y = 0;
    iter_count = 0;

    last_pos_x1 = 0;
    last_pos_y1 = 0;
    last_theta1 = 0;
    
    clf();
end


subplot(2,1,1);
dt = 1.0/settings.sensorPollFreq;
[ disp_x, disp_y ] = calculate_fly_position_no_yaw(vel_forward, vel_side, dt, last_pos_x, last_pos_y);
hold on;
plot(disp_x, disp_y, 'color', 'blue');
xlabel('X displacement (au)');
ylabel('Y displacement (au)');
title('No yaw');
last_pos_x = disp_x(end);
last_pos_y = disp_y(end);

subplot(2,1,2);
[disp_x, disp_y, theta] = calculate_fly_position_with_yaw(vel_forward, vel_side, vel_yaw, dt, last_pos_x1, last_pos_y1, last_theta1);
hold on;
plot(disp_x, disp_y, 'color', 'blue');
xlabel('X displacement (au)');
ylabel('Y displacement (au)');
title('With yaw');

last_pos_x1 = disp_x(end);
last_pos_y1 = disp_y(end);
last_theta1 = theta(end);
iter_count = iter_count + 1;

end

