function queueMoreData_2Piezos(src, event, fid, chunklength)

outData = fread(fid, [2, chunklength*8] ,'single' ); % 8 bits for each timepoint
if ~feof(fid) 
    queueOutputData(src, outData');
elseif length(outData) > chunklength/3
    queueOutputData(src, outData');
else % will discard the last bits of data
    disp('experiment all done')
    src.stop;
end

end