function varargout = fly_tracker_TO(varargin)
% FLY_TRACKER_TO MATLAB code for fly_tracker_TO.fig
%      FLY_TRACKER_TO, by itself, creates a new FLY_TRACKER_TO or raises the existing
%      singleton*.
%
%      H = FLY_TRACKER_TO returns the handle to a new FLY_TRACKER_TO or the handle to
%      the existing singleton*.
%
%      FLY_TRACKER_TO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLY_TRACKER_TO.M with the given input arguments.
%
%      FLY_TRACKER_TO('Property','Value',...) creates a new FLY_TRACKER_TO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fly_tracker_TO_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fly_tracker_TO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fly_tracker_TO

% Last Modified by GUIDE v2.5 08-Aug-2017 08:55:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @fly_tracker_TO_OpeningFcn, ...
    'gui_OutputFcn',  @fly_tracker_TO_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before fly_tracker_TO is made visible.
function fly_tracker_TO_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fly_tracker_TO (see VARARGIN)

% Choose default command line output for fly_tracker_TO
handles.output = hObject;

% UIWAIT makes fly_tracker_TO wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Prompt for an experiment directory
dname = uigetdir('D:\Data\', 'Please chose an experiment directory.');
handles.experiment_dir = dname;

ghandles = guihandles(hObject);
set(ghandles.experiment_dir_edit, 'String', dname);

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = fly_tracker_TO_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in experiment_dir_button.
function experiment_dir_button_Callback(hObject, eventdata, handles)
% hObject    handle to experiment_dir_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dname = uigetdir('C:\Users\sasha\Dropbox\Wilson_lab\paper_1\data\descending_neurons\');
handles.experiment_dir = dname;

ghandles = guihandles(hObject);
set(ghandles.experiment_dir_edit, 'String', dname);

% Update handles structure
guidata(hObject, handles);



function experiment_dir_edit_Callback(hObject, eventdata, handles)
% hObject    handle to experiment_dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of experiment_dir_edit as text
%        str2double(get(hObject,'String')) returns contents of experiment_dir_edit as a double


% --- Executes during object creation, after setting all properties.
function experiment_dir_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to experiment_dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pre_stim_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pre_stim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pre_stim_edit as text
%        str2double(get(hObject,'String')) returns contents of pre_stim_edit as a double


% --- Executes during object creation, after setting all properties.
function pre_stim_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pre_stim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stim_edit_Callback(hObject, eventdata, handles)
% hObject    handle to stim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim_edit as text
%        str2double(get(hObject,'String')) returns contents of stim_edit as a double


% --- Executes during object creation, after setting all properties.
function stim_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function post_stim_edit_Callback(hObject, eventdata, handles)
% hObject    handle to post_stim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of post_stim_edit as text
%        str2double(get(hObject,'String')) returns contents of post_stim_edit as a double


% --- Executes during object creation, after setting all properties.
function post_stim_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to post_stim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function inter_trial_period_edit_Callback(hObject, eventdata, handles)
% hObject    handle to inter_trial_period_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inter_trial_period_edit as text
%        str2double(get(hObject,'String')) returns contents of inter_trial_period_edit as a double


% --- Executes during object creation, after setting all properties.
function inter_trial_period_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inter_trial_period_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in stim_type_menu.
function stim_type_menu_Callback(hObject, eventdata, handles)
% hObject    handle to stim_type_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns stim_type_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stim_type_menu


% --- Executes during object creation, after setting all properties.
function stim_type_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_type_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in using_2p_button.
function using_2p_button_Callback(hObject, eventdata, handles)
% hObject    handle to using_2p_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of using_2p_button



function task_file_edit_Callback(hObject, eventdata, handles)
% hObject    handle to task_file_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of task_file_edit as text
%        str2double(get(hObject,'String')) returns contents of task_file_edit as a double


% --- Executes during object creation, after setting all properties.
function task_file_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to task_file_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in task_file_button.
function task_file_button_Callback(hObject, eventdata, handles)
% hObject    handle to task_file_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uigetfile('D:\TaskFiles\*.txt','Select a task file');

handles.taskfile_path = [PathName '\' FileName];

ghandles = guihandles(hObject);
set(ghandles.task_file_edit, 'String', handles.taskfile_path);

guidata(hObject, handles);

% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ghandles = guihandles(hObject);

task_filepath = get(ghandles.task_file_edit, 'String');
if( isempty(task_filepath) == 1 )
    errordlg('Please set the task file path.')
    return;
end

run_obj.pre_stim_t = str2num(get(ghandles.pre_stim_edit, 'String'));
run_obj.stim_t = str2num(get(ghandles.stim_edit, 'String'));
run_obj.post_stim_t = str2num(get(ghandles.post_stim_edit, 'String'));
run_obj.inter_trial_t = str2num(get(ghandles.inter_trial_period_edit, 'String'));

contents = get(ghandles.stim_type_menu,'String');
run_obj.stim_type = contents; % contents{get(ghandles.stim_type_menu,'Value')}; no longer a dropdown menu

run_obj.experiment_dir = handles.experiment_dir;
run_obj.session_id = str2num(get(ghandles.session_id_edit, 'String'));
run_obj.sessiod_id_hdl = ghandles.session_id_edit;
run_obj.taskfile_path = task_filepath;

run_obj.injection_current = str2num(get(ghandles.external_command_edit, 'String'));

run_obj.patch_id = str2num(get(ghandles.patch_id_edit, 'String'));

run_obj.IsBehavior = get(ghandles.checkbox_behavior,'Value');
%start_trials(run_obj);
start_trials_continuous_TO(run_obj);

guidata(hObject, handles);

function session_id_edit_Callback(hObject, eventdata, handles)
% hObject    handle to session_id_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of session_id_edit as text
%        str2double(get(hObject,'String')) returns contents of session_id_edit as a double


% --- Executes during object creation, after setting all properties.
function session_id_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to session_id_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function external_command_edit_Callback(hObject, eventdata, handles)
% hObject    handle to external_command_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of external_command_edit as text
%        str2double(get(hObject,'String')) returns contents of external_command_edit as a double


% --- Executes during object creation, after setting all properties.
function external_command_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to external_command_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stop_button.
function stop_button_Callback(hObject, eventdata, handles)
% hObject    handle to stop_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global session_obj;
session_obj.stop();
release( session_obj );
disp(['Stopped acquisition.']);



function patch_id_edit_Callback(hObject, eventdata, handles)
% hObject    handle to patch_id_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of patch_id_edit as text
%        str2double(get(hObject,'String')) returns contents of patch_id_edit as a double


% --- Executes during object creation, after setting all properties.
function patch_id_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to patch_id_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in record_pipette_resistance_button.
function record_pipette_resistance_button_Callback(hObject, eventdata, handles)
% hObject    handle to record_pipette_resistance_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ghandles = guihandles(hObject);

Ch = rig_params_TO;
s = daq.createSession('ni');
ai_channels_used = [0:15];
aI = s.addAnalogInputChannel(Ch.DevID, ai_channels_used, 'Voltage');
for i=1:length(ai_channels_used)
    aI(i).InputType = 'SingleEnded';
end

settings = sensor_settings;

SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
s.DurationInSeconds = 1.0;

[trial_data, trial_time] = s.startForeground();
release(s);

[ current, voltage, scaled ] = get_scaled_voltage_and_current_TO( trial_data );

pipette_resistance = pipetteResistanceCalc( current);

set(ghandles.pipette_resistance_text, 'String', [num2str(pipette_resistance, 3) ' MOhm']);

patch_id = str2num(get(ghandles.patch_id_edit, 'String'));

experiment_dir = handles.experiment_dir;

pr_filename = [experiment_dir '/pipette_resistance_patch_A_' num2str(patch_id) '.txt'];


fileID = fopen(pr_filename,'w+');
fprintf(fileID,'%f MOhm\n',pipette_resistance);
fclose(fileID);

pr_filename_mat = [experiment_dir '/pipette_resistance_patch_A_' num2str(patch_id) '.mat'];
save(pr_filename_mat, 'current', 'trial_time');

guidata(hObject, handles);

% --- Executes on button press in record_access_resistance_button.
function record_access_resistance_button_Callback(hObject, eventdata, handles)
% hObject    handle to record_access_resistance_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ghandles = guihandles(hObject);

Ch = rig_params_TO;
s = daq.createSession('ni');
ai_channels_used = [0:15];
aI = s.addAnalogInputChannel(Ch.DevID, ai_channels_used, 'Voltage');
for i=1:length(ai_channels_used)
    aI(i).InputType = 'SingleEnded';
end

settings = sensor_settings;

settings.sampRate = 40000; % higher sampling rate for getting the peak
SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
s.DurationInSeconds = 1.0;

[trial_data, trial_time] = s.startForeground();
release(s);

[ current, voltage, scaled ] = get_scaled_voltage_and_current_TO( trial_data );

%access_resistance = accessResistanceCalc( current, SAMPLING_RATE );
[R_access, R_membrane] = accessResistanceCalc_TO( current, voltage, trial_time, SAMPLING_RATE );

set(ghandles.access_resistance_text, 'String', ['R_a: ', num2str(R_access,3) ' MOhm']);
set(ghandles.membrane_resistance_text, 'String', ['R_in: ', num2str(R_membrane,3) ' MOhm']);

patch_id = str2num(get(ghandles.patch_id_edit, 'String'));

experiment_dir = handles.experiment_dir;

pr_filename = [experiment_dir '/access_resistance_patch_A_' num2str(patch_id) '.txt'];

fileID = fopen(pr_filename,'w+');
fprintf(fileID,'%f MOhm\n',R_access);
fclose(fileID);

pr_filename_mat = [experiment_dir '/access_resistance_patch_A_' num2str(patch_id) '.mat'];
save(pr_filename_mat, 'current', 'trial_time');

guidata(hObject, handles);

% --- Executes on button press in record_i_0_resting_voltage.
function record_i_0_resting_voltage_Callback(hObject, eventdata, handles)
% hObject    handle to record_i_0_resting_voltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ghandles = guihandles(hObject);

Ch = rig_params_TO;
s = daq.createSession('ni');
ai_channels_used = [0:15];
aI = s.addAnalogInputChannel(Ch.DevID, ai_channels_used, 'Voltage');
for i=1:length(ai_channels_used)
    aI(i).InputType = 'SingleEnded';
end

settings = sensor_settings;

SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
s.DurationInSeconds = 1.0;

[trial_data, trial_time] = s.startForeground();
release(s);

[ current, voltage, scaled ] = get_scaled_voltage_and_current_TO( trial_data );

resting_voltage = mean( voltage );

set(ghandles.i_0_resting_voltage_text, 'String', [num2str(resting_voltage,3) ' mV']);

patch_id = str2num(get(ghandles.patch_id_edit, 'String'));

experiment_dir = handles.experiment_dir;

pr_filename = [experiment_dir '/initial_resting_voltage_patch_A_' num2str(patch_id) '.txt'];

fileID = fopen(pr_filename,'w+');
fprintf(fileID, '%f mV\n', resting_voltage);
fclose(fileID);

pr_filename_mat = [experiment_dir '/initial_resting_voltage_patch_A_' num2str(patch_id) '.mat'];
save(pr_filename_mat, 'voltage', 'trial_time');

guidata(hObject, handles);


% --- Executes on button press in init_ball_button.
function init_ball_button_Callback(hObject, eventdata, handles)
% hObject    handle to init_ball_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This requires the ball to not be moving.
init_ball_parameters_TO(handles.experiment_dir);


% --- Executes on button press in record_pipette_resistance_button_B.
function record_pipette_resistance_button_B_Callback(hObject, eventdata, handles)
% hObject    handle to record_pipette_resistance_button_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ghandles = guihandles(hObject);

Ch = rig_params_TO;
s = daq.createSession('ni');
ai_channels_used = [0:15];
aI = s.addAnalogInputChannel(Ch.DevID, ai_channels_used, 'Voltage');
for i=1:length(ai_channels_used)
    aI(i).InputType = 'SingleEnded';
end

settings = sensor_settings;

SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
s.DurationInSeconds = 1.0;

[trial_data, trial_time] = s.startForeground();
release(s);

[ current, voltage, scaled ] = get_scaled_voltage_and_current_TO( trial_data );

% figure;
% plot(trial_time, current);

pipette_resistance = pipetteResistanceCalc( current );

set(ghandles.pipette_resistance_text_B, 'String', [num2str(pipette_resistance) ' MOhm']);

patch_id = str2num(get(ghandles.patch_id_edit_B, 'String'));

experiment_dir = handles.experiment_dir;

pr_filename = [experiment_dir '/pipette_resistance_patch_B_' num2str(patch_id) '.txt'];

fileID = fopen(pr_filename,'w+');
fprintf(fileID,'%f MOhm\n',pipette_resistance);
fclose(fileID);

pr_filename_mat = [experiment_dir '/pipette_resistance_patch_B_' num2str(patch_id) '.mat'];
save(pr_filename_mat, 'current', 'trial_time');

guidata(hObject, handles);

% --- Executes on button press in record_access_resistance_button_B.
function record_access_resistance_button_B_Callback(hObject, eventdata, handles)
% hObject    handle to record_access_resistance_button_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ghandles = guihandles(hObject);

Ch = rig_params_TO;
s = daq.createSession('ni');
aI = s.addAnalogInputChannel(Ch.DevID, [0:15], 'Voltage');
for i=1:length(ai_channels_used)
    aI(i).InputType = 'SingleEnded';
end

settings = sensor_settings;

s.Rate = 40000; % higher sampling rate to capture the transients
s.DurationInSeconds = 1.0;

[trial_data, trial_time] = s.startForeground();
release(s);

[ current, voltage, scaled ] = get_scaled_voltage_and_current_TO( trial_data );

access_resistance = accessResistanceCalc( current, SAMPLING_RATE );

set(ghandles.access_resistance_text_B, 'String', [num2str(access_resistance) ' MOhm']);

patch_id = str2num(get(ghandles.patch_id_edit_B, 'String'));

experiment_dir = handles.experiment_dir;

pr_filename = [experiment_dir '/access_resistance_patch_B_' num2str(patch_id) '.txt'];

fileID = fopen(pr_filename,'w+');
fprintf(fileID,'%f MOhm\n',access_resistance);
fclose(fileID);

pr_filename_mat = [experiment_dir '/access_resistance_patch_B_' num2str(patch_id) '.mat'];
save(pr_filename_mat, 'current', 'trial_time');

guidata(hObject, handles);

% --- Executes on button press in record_i_0_resting_voltage_B.
function record_i_0_resting_voltage_B_Callback(hObject, eventdata, handles)
% hObject    handle to record_i_0_resting_voltage_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ghandles = guihandles(hObject);

Ch = rig_params_TO;
s = daq.createSession('ni');
aI = s.addAnalogInputChannel(Ch.DevID, [0:15], 'Voltage');

settings = sensor_settings;

SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
s.DurationInSeconds = 1.0;

[trial_data, trial_time] = s.startForeground();
release(s);

[ current, voltage, scaled ] = get_scaled_voltage_and_current_TO( trial_data );

resting_voltage = mean( voltage );

set(ghandles.i_0_resting_voltage_text_B, 'String', [num2str(resting_voltage) ' mV']);

patch_id = str2num(get(ghandles.patch_id_edit_B, 'String'));

experiment_dir = handles.experiment_dir;

pr_filename = [experiment_dir '/initial_resting_voltage_patch_B_' num2str(patch_id) '.txt'];

fileID = fopen(pr_filename,'w+');
fprintf(fileID, '%f mV\n', resting_voltage);
fclose(fileID);

pr_filename_mat = [experiment_dir '/initial_resting_voltage_patch_B_' num2str(patch_id) '.mat'];
save(pr_filename_mat, 'voltage', 'trial_time');

guidata(hObject, handles);

function patch_id_edit_B_Callback(hObject, eventdata, handles)
% hObject    handle to patch_id_edit_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of patch_id_edit_B as text
%        str2double(get(hObject,'String')) returns contents of patch_id_edit_B as a double


% --- Executes during object creation, after setting all properties.
function patch_id_edit_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to patch_id_edit_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in button_group_behavior.
function button_group_behavior_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in button_group_behavior
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ghandles = guihandles(hObject);
if ghandles.radiobutton_OGA.Value==1
    set(ghandles.pre_stim_edit','String',num2str(4))
    set(ghandles.pre_stim_edit,'Enable','off');
    set(ghandles.stim_edit','String',num2str(0.5))
    set(ghandles.stim_edit,'Enable','off');
    set(ghandles.post_stim_edit','String',num2str(4));
    set(ghandles.post_stim_edit,'Enable','off');
    set(ghandles.inter_trial_period_edit','String',num2str(3));
    set(ghandles.inter_trial_period_edit','Enable','off');
elseif ghandles.radiobutton_OTT.Value==1
    set(ghandles.pre_stim_edit,'Enable','on');
    set(ghandles.stim_edit,'Enable','on');
    set(ghandles.post_stim_edit,'Enable','on');
    set(ghandles.inter_trial_period_edit','Enable','on');
    set(ghandles.pre_stim_edit','String',num2str(3));
    set(ghandles.stim_edit','String',num2str(0.5))
    set(ghandles.post_stim_edit','String',num2str(3));
    set(ghandles.inter_trial_period_edit','String',num2str(5));
end


% --- Executes on button press in record_activity.
function record_activity_Callback(hObject, eventdata, handles)
% hObject    handle to record_activity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ghandles = guihandles(hObject);
set(ghandles.record_activity,'ForegroundColor','r')

Ch = rig_params_TO;
s = daq.createSession('ni');
ai_channels_used = [0:15];
aI = s.addAnalogInputChannel(Ch.DevID, ai_channels_used, 'Voltage');
for i=1:length(ai_channels_used)
    aI(i).InputType = 'SingleEnded';
end

settings = sensor_settings;

s.Rate = 40000; % higher sampling rate to capture transients
s.DurationInSeconds = str2num(get(ghandles.edit_duration,'String'));

[trial_data, trial_time] = s.startForeground();
release(s);

set(ghandles.record_activity,'ForegroundColor','k')

[ current, voltage, scaled, mode] = get_scaled_voltage_and_current_TO( trial_data );
%[ current, voltage] = get_scaled_voltage_and_current_A( trial_data );

figure(20); clf;
s1=subplot(311);
plot(trial_time,scaled,'r')
set(gca,'color','none')
box off
if strcmp(mode,'I')
    ylabel('Current (pA)','fontsize',12)
elseif strcmp(mode,'V')
    ylabel('Voltage (mV)','fontsize',12)
end
s2=subplot(312);
plot(trial_time,current,'b')
ylabel('Current (pA)','fontsize',12)
box off
set(gca,'color','none')
s3=subplot(313);
plot(trial_time,voltage,'k')
box off
set(gca,'color','none')
ylabel('Voltage (mV)','fontsize',12)
xlabel('Time (s)','fontsize',12)
linkaxes([s1,s2],'x')

experiment_dir = handles.experiment_dir;

pr_filename_mat = [experiment_dir '/activity_', datestr(now,'yyyymmdd_HHMMSS'), '.mat'];
save(pr_filename_mat, 'current', 'voltage', 'trial_time','trial_data');



function edit_duration_Callback(hObject, eventdata, handles)
% hObject    handle to edit_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_duration as text
%        str2double(get(hObject,'String')) returns contents of edit_duration as a double


% --- Executes during object creation, after setting all properties.
function edit_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in record_seal_resistance_button.
function record_seal_resistance_button_Callback(hObject, eventdata, handles)
% hObject    handle to record_seal_resistance_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ghandles = guihandles(hObject);

Ch = rig_params_TO;
s = daq.createSession('ni');
ai_channels_used = [0:15];
aI = s.addAnalogInputChannel(Ch.DevID, ai_channels_used, 'Voltage');
for i=1:length(ai_channels_used)
    aI(i).InputType = 'SingleEnded';
end

settings = sensor_settings;

settings.sampRate = 40000; % higher sampling rate for getting the peak
SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
s.DurationInSeconds = 1.0; %%%%%%%%%%%

set(handles.record_seal_resistance_button,'ForegroundColor',[1,0,0]) % change the font to red
[trial_data, trial_time] = s.startForeground();
release(s);
set(handles.record_seal_resistance_button,'ForegroundColor',[0,0,0]) % change the font to red

 [ current, voltage ] = get_scaled_voltage_and_current_A( trial_data );
%[ current, voltage ] = get_nonscaled_voltage_and_current_TO( trial_data );
Rseal = sealResistanceCalc_TO(current, voltage, trial_time);

set(ghandles.seal_resistance_text, 'String', [num2str(Rseal/1e9, 3) ' GOhm']);

patch_id = str2num(get(ghandles.patch_id_edit, 'String'));

experiment_dir = handles.experiment_dir;

pr_filename = [experiment_dir '/seal_resistance_patch_A_' num2str(patch_id) '.txt'];

fileID = fopen(pr_filename,'w+');
fprintf(fileID,'%f GOhm\n',Rseal/1e9);
fclose(fileID);

pr_filename_mat = [experiment_dir '/seal_resistance_patch_A_' num2str(patch_id) '.mat'];
save(pr_filename_mat, 'current','voltage', 'trial_time');

guidata(hObject, handles);



% --- Executes on button press in checkbox_behavior.
function checkbox_behavior_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_behavior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_behavior
