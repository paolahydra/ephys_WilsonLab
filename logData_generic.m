function logData_generic(src, evt, fid)
% Add the time stamp and the data values to data. To write data sequentially,
% transpose the matrix.

%   Copyright 2011 The MathWorks, Inc.

data = [evt.TimeStamps, evt.Data]' ;
fwrite(fid,data,'double');

% puoi dare come input degli assi di una finestra che hai aperto prima in
% cui plottare dati online ogni volta che sono disponibili (hold off,
% oppure usando le timestamps)

end