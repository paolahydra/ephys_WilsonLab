function [] = display_trial( task, trial_time, trial_data, viz_figs )

colors = {'red', 'green', 'blue', 'black'};
cur_color = '';

if( (strcmp(task, 'LeftOdor') == 1 ) )
    cur_color = colors{1};
elseif( strcmp(task, 'RightOdor') == 1 )
    cur_color = colors{2};
elseif( strcmp(task, 'BothOdor') == 1 )
    cur_color = colors{3};
elseif( strcmp(task, 'NaturalOdor') == 1 )
    cur_color = colors{4};
else
    disp(['ERROR: Task: ' task ' is not recognized.']);
end

settings = sensor_settings;

[ t, vel_forward, vel_side, vel_yaw ] = get_velocity_ephysrig(trial_time, trial_data); 

% Display trial velocities
%figure(viz_figs.velocity_tc_fig);
set( 0, 'CurrentFigure', viz_figs.velocity_tc_fig )

% Plot forward
subplot(3,1,1);
hold on;
plot( t, vel_forward, 'color', cur_color );
ylabel('au/s');
xlabel('Time (s)');
xlim([0 trial_time(end)]);
title('Forward velocity');

subplot(3,1,2);
hold on;
plot( t, vel_side, 'color', cur_color );
ylabel('au/s');
xlabel('Time (s)');
xlim([0 trial_time(end)]);
title('Side velocity');

subplot(3,1,3);
hold on;
plot( t, vel_yaw, 'color', cur_color );
ylabel('au/s');
xlabel('Time (s)');
xlim([0 trial_time(end)]);
title('Yaw velocity');

% Display trial raw trajectory
% figure(viz_figs.run_traj_fig);
set( 0, 'CurrentFigure', viz_figs.run_traj_fig )

dt = 1.0/settings.sensorPollFreq;
%[disp_x, disp_y, theta] = calculate_fly_position_with_yaw(vel_forward, vel_side, vel_yaw, dt, 0, 0, 0);
[disp_x, disp_y] = calculate_fly_position_no_yaw(vel_forward, vel_side, dt, 0, 0);
hold on;
plot(disp_x, disp_y, 'color', cur_color);
xlabel('X displacement (au)');
ylabel('Y displacement (au)');

%h = legend('Both Odor', 'Left Odor', 'Right Odor');
%legendlinecolors(h, {'blue', 'red', 'green'});
h = legend('Left Odor', 'Right Odor');
legendlinecolors(h, {'red', 'green'});

%figure( viz_figs.ephys_fig );
set( 0, 'CurrentFigure', viz_figs.ephys_fig )

current = trial_data(:,1);
voltage = trial_data(:,2);

subplot(2,1,1);
hold on;
plot(trial_time, current);
title( 'Current' );

subplot(2,1,2);
VOLTAGE_SCALING_FACTOR = 10;
hold on;
plot(trial_time, VOLTAGE_SCALING_FACTOR*voltage);
title( 'Voltage' );

end

