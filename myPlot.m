function[h]= myPlot(Data,varargin)

if length(Data) < 3
    disp('invalid Print, length<3');
    return;
end

[col, row] = size(Data);
if isempty(intersect([1,2,3,4,6],row))
    Data = Data.';
end
[col, row] = size(Data);

if 1 == row
    h = plot(Data, varargin{:});
elseif 2 == row
    h = plot(Data(:,1), Data(:,2), varargin{:});
elseif 3 == row
    h = plot3(Data(:,1), Data(:,2), Data(:,3), varargin{:});
elseif 4 == row
    h = quiver(Data(:,1), Data(:,2), Data(:,3), Data(:,4), varargin{:});
elseif 6 == row
    h = quiver3(Data(:,1), Data(:,2), Data(:,3), Data(:,4), Data(:,5), Data(:,6), varargin{:});
else
    fprintf('invalid Print! col=%d row=%d\n', col, row);
end


