function made = MakeDirIfMissing(dirPath)

% MakeDirIfMissing create a directory if it doesnot exist. The underlying
% call is matlab built-in function mkdir.
% dirPath: string, full path of a directory to create
% Returned made is a logical scalar. true means the specified dirPath is
% created in this call, while false means the dirPath has already been
% there.
%
% By Ji Zhou

made = false;
if ~exist(dirPath, 'dir')
    mkdir(dirPath);
    made = true;
end

end