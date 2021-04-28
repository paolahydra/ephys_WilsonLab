function varargout = fly_tracker(varargin)
% FLY_TRACKER MATLAB code for fly_tracker.fig
%      FLY_TRACKER, by itself, creates a new FLY_TRACKER or raises the existing
%      singleton*.
%
%      H = FLY_TRACKER returns the handle to a new FLY_TRACKER or the handle to
%      the existing singleton*.
%
%      FLY_TRACKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLY_TRACKER.M with the given input arguments.
%
%      FLY_TRACKER('Property','Value',...) creates a new FLY_TRACKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fly_tracker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fly_tracker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fly_tracker

% Last Modified by GUIDE v2.5 30-Nov-2016 09:58:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fly_tracker_OpeningFcn, ...
                   'gui_OutputFcn',  @fly_tracker_OutputFcn, ...
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


% --- Executes just before fly_tracker is made visible.
function fly_tracker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fly_tracker (see VARARGIN)

% Choose default command line output for fly_tracker
handles.output = hObject;

% UIWAIT makes fly_tracker wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Prompt for an experiment directory
dname = uigetdir('C:\Users\sasha\Dropbox\Wilson_lab\paper_1\data\descending_neurons\', 'Please chose an experiment directory.');
handles.experiment_dir = dname;

ghandles = guihandles(hObject);
set(ghandles.experiment_dir_edit, 'String', dname);

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = fly_tracker_OutputFcn(hObject, eventdata, handles) 
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

[FileName,PathName] = uigetfile('C:\Users\sasha\Dropbox\Wilson_lab\fly_yaw_reader_ephys\task_files\*.txt','Select a task file');

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
run_obj.stim_type = contents; % contents{get(ghandles.stim_type_menu,'Value')};

run_obj.experiment_dir = handles.experiment_dir;
run_obj.session_id = str2num(get(ghandles.session_id_edit, 'String'));
run_obj.sessiod_id_hdl = ghandles.session_id_edit;
run_obj.taskfile_path = task_filepath;

run_obj.injection_current = str2num(get(ghandles.external_command_edit, 'String'));

run_obj.patch_id = str2num(get(ghandles.patch_id_edit, 'String'));

%start_trials(run_obj);
start_trials_continuous(run_obj);

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

s = daq.createSession('ni');
aI = s.addAnalogInputChannel('Dev3', [0:31], 'Voltage');

settings = sensor_settings;

SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
s.DurationInSeconds = 1.0;

[trial_data, trial_time] = s.startForeground();
release(s);

[ current, voltage ] = get_scaled_voltage_and_current_A( trial_data );

pipette_resistance = pipetteResistanceCalc( current );

set(ghandles.pipette_resistance_text, 'String', [num2str(pipette_resistance) ' MOhm']); 

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

s = daq.createSession('ni');
aI = s.addAnalogInputChannel('Dev3', [0:31], 'Voltage');

settings = sensor_settings;

SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
s.DurationInSeconds = 1.0;

[trial_data, trial_time] = s.startForeground();
release(s);

[ current, voltage ] = get_scaled_voltage_and_current_A( trial_data );

access_resistance = accessResistanceCalc( current, SAMPLING_RATE );

set(ghandles.access_resistance_text, 'String', [num2str(access_resistance) ' MOhm']); 

patch_id = str2num(get(ghandles.patch_id_edit, 'String'));

experiment_dir = handles.experiment_dir;

pr_filename = [experiment_dir '/access_resistance_patch_A_' num2str(patch_id) '.txt'];

fileID = fopen(pr_filename,'w+');
fprintf(fileID,'%f MOhm\n',access_resistance);
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

s = daq.createSession('ni');
aI = s.addAnalogInputChannel('Dev3', [0:31], 'Voltage');

settings = sensor_settings;

SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
s.DurationInSeconds = 1.0;

[trial_data, trial_time] = s.startForeground();
release(s);

[ current, voltage ] = get_scaled_voltage_and_current_A( trial_data );

resting_voltage = mean( voltage );

set(ghandles.i_0_resting_voltage_text, 'String', [num2str(resting_voltage) ' mV']); 

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
init_ball_parameters(handles.experiment_dir);


% --- Executes on button press in record_pipette_resistance_button_B.
function record_pipette_resistance_button_B_Callback(hObject, eventdata, handles)
% hObject    handle to record_pipette_resistance_button_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ghandles = guihandles(hObject);

s = daq.createSession('ni');
ai_channels_used = [0:31];
aI = s.addAnalogInputChannel('Dev3', ai_channels_used, 'Voltage');
for i=1:length(ai_channels_used)
    aI(i).InputType = 'SingleEnded';
end

settings = sensor_settings;

SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
s.DurationInSeconds = 1.0;

[trial_data, trial_time] = s.startForeground();
release(s);

[ current, voltage ] = get_scaled_voltage_and_current_B( trial_data );

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

s = daq.createSession('ni');
aI = s.addAnalogInputChannel('Dev3', [0:31], 'Voltage');

settings = sensor_settings;

SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
s.DurationInSeconds = 1.0;

[trial_data, trial_time] = s.startForeground();
release(s);

[ current, voltage ] = get_scaled_voltage_and_current_B( trial_data );

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

s = daq.createSession('ni');
aI = s.addAnalogInputChannel('Dev3', [0:31], 'Voltage');

settings = sensor_settings;

SAMPLING_RATE = settings.sampRate;
s.Rate = SAMPLING_RATE;
s.DurationInSeconds = 1.0;

[trial_data, trial_time] = s.startForeground();
release(s);

[ current, voltage ] = get_scaled_voltage_and_current_B( trial_data );

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
