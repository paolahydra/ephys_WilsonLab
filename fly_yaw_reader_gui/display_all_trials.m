function [] = display_all_trials( task, trial_time, trial_data, viz_figs, pre_stim_t, stim_t, experiment_dir )

colors = { rgb('Red'), rgb('Green'), rgb('Blue'), rgb('Black'), rgb('Brown'), rgb('Purple')};
cur_color = '';

if( (strcmp(task, 'LeftOdor') == 1 ) || (strcmp(task, 'WideFieldLight') == 1 ) )
    cur_color = colors{1};
elseif( strcmp(task, 'RightOdor') == 1 )
    cur_color = colors{2};
elseif( strcmp(task, 'BothOdor') == 1 )
    cur_color = colors{3};
elseif( strcmp(task, 'NaturalOdor') == 1 )
    cur_color = colors{4};
elseif( strcmp(task, 'ExternalCommandDepol') == 1 )
    cur_color = colors{5};
elseif( strcmp(task, 'ExternalCommandHypopol') == 1 )
    cur_color = colors{6};
elseif( strcmp(task, 'OdorNoWind') == 1 )
    cur_color = colors{4};
elseif( strcmp(task, 'OdorLeftWind') == 1 )
    cur_color = colors{1};
elseif( strcmp(task, 'OdorRightWind') == 1 )
    cur_color = colors{2};
else
    disp(['ERROR: Task: ' task ' is not recognized.']);
end

use_calibration = 1;
[ t, vel_forward, vel_side, vel_yaw ] = get_velocity_ephysrig(trial_time, trial_data, experiment_dir, use_calibration); 

% Display trial raw trajectory
% figure(viz_figs.run_traj_fig);
set( 0, 'CurrentFigure', viz_figs.run_traj_fig )

settings = sensor_settings;
dt = 1.0/settings.sensorPollFreq;
%[disp_x, disp_y, theta] = calculate_fly_position_with_yaw(vel_forward, vel_side, vel_yaw, dt, 0, 0, 0);
[disp_x, disp_y] = calculate_fly_position_no_yaw(vel_forward, vel_side, dt, 0, 0);
hold on;
plot(disp_x, disp_y, 'color', cur_color);
xlabel('X displacement (mm)');
ylabel('Y displacement (mm)');

set( 0, 'CurrentFigure', viz_figs.all_trials_fig )

% Plot forward
subplot(4,1,1);
hold on;
plot( t, vel_forward, 'color', cur_color );
ylabel('Fwd velocity (au/s)');
xlim([0 trial_time(end)]);

t_vel_first = t(1);

yy = ylim;
y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
hh = fill([ t_vel_first+pre_stim_t t_vel_first+pre_stim_t (t_vel_first+pre_stim_t+stim_t) (t_vel_first+pre_stim_t+stim_t) ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

subplot(4,1,2);
hold on;
plot( t, vel_yaw, 'color', cur_color );
ylabel('Yaw velocity (au/s)');
xlim([0 trial_time(end)]);
yy = ylim;
y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
hh = fill([ t_vel_first+pre_stim_t t_vel_first+pre_stim_t (t_vel_first+pre_stim_t+stim_t) (t_vel_first+pre_stim_t+stim_t) ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

[currentA, voltageA, currentB, voltageB] = get_dual_scaled_voltage_and_current( trial_data );

t_first = trial_time(1);

subplot(4,1,3);
hold on;
plot(trial_time, currentA, 'color', cur_color );
plot(trial_time, currentB, 'color', cur_color, 'LineStyle', '--' );
xlim([0 trial_time(end)]);
ylabel('Current (pA)');
yy = ylim;
y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
hh = fill([ t_first+pre_stim_t t_first+pre_stim_t (t_first+pre_stim_t+stim_t) (t_first+pre_stim_t+stim_t) ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

subplot(4,1,4);
hold on;
plot(trial_time, voltageA, 'color', cur_color );
plot(trial_time, voltageB, 'color', cur_color, 'LineStyle', '--' );
xlim([0 trial_time(end)]);
xlabel('Time (s)');
ylabel('Voltage (mV)')
yy = ylim;
y_min = yy(1)-yy(1)*0.01; y_max = yy(2);
hh = fill([ t_first+pre_stim_t t_first+pre_stim_t (t_first+pre_stim_t+stim_t) (t_first+pre_stim_t+stim_t) ],[y_min y_max y_max y_min ], rgb('Wheat'));
set(gca,'children',circshift(get(gca,'children'),-1));
set(hh, 'EdgeColor', 'None');

end

