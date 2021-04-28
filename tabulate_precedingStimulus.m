function counts = tabulate_precedingStimulus(x)
%TABULATE Frequency table.
%   TABLE = TABULATE(X) takes a vector X and returns a square matrix, TABLE.
%   For each value in X (rows), it return the count of every value
%   (columns) found in the position preceding it.

   
%   Copyright 1993-2011 The MathWorks, Inc.


isnum = isnumeric(x);
if isnum && ~isfloat(x)
    % use of hist() below requires float
    x = double(x);
end
if isnum
   if min(size(x)) > 1,
      error(message('stats:tabulate:InvalidData'));
   end

   y = x(~isnan(x));
else
   y = x;
end

if ~isnum || any(y ~= round(y)) || any(y < 1)
   docell = true;
   [y,yn,yl] = grp2idx(y);
   maxlevels = length(yn);
else
   docell = false;
   maxlevels = max(y);
   %yn = cellstr(num2str((1:maxlevels)'));
end

[~, values] = hist(y,(1:maxlevels));

counts = zeros(length(values));
y = y(:);
slidedY = [nan; y(1:end-1)];
for c = 1:length(values)
    precvalues = slidedY(y == values(c));
    for i = 1:length(values)
        counts(c,i) = sum(precvalues==values(i));
    end
end



end