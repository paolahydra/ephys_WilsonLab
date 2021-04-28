function [ tasks ] = read_task_file( task_file_path )

tasks = {};

fid = fopen(task_file_path);
tline = fgetl(fid);

i=1;
while ischar(tline)
    tasks{i} = tline;
    %disp(tline);
    tline = fgetl(fid);
    i = i + 1;
end

fclose(fid);

end

