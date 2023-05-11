function results = ReadYaml(filePath)
% Lloyd Russell 2017
% Simple little function to read simple little YAML format parameter files

% read file line by line
fid = fopen(filePath, 'r');
data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
fclose(fid);

% remove empty lines
data = deblank(data{1});
data(cellfun('isempty', data)) = [];

% prepare final results structure
results = [];

catData = 0;
setProp = 0;
% parse the contents (line by line)

for i = 1:numel(data)
    
    % extract this line
    thisLine = data{i};
    
    % ignore if this line is a comment
    if strcmpi(thisLine(1), '#') ||strcmpi(thisLine(1), '%') || strcmpi(thisLine(1), '-')
        continue
    end
    
    % find the seperator between key and value
    sepIndex = find(thisLine==':', 1, 'first');
    
    % get the key name (remove whitespace)
    key = strtrim(thisLine(1:sepIndex-1));  
    
    if isempty(key) && catData
        results.(lastKey1).(lastKey2) = [results.(lastKey1).(lastKey2) thisLine];
        if strcmpi(thisLine(end),']')
            catData = 0;
            results.(lastKey1).(lastKey2) = str2double(strsplit(results.(lastKey1).(lastKey2)(2:end-1),','));
            if isfield(results.(lastKey1),'rows') 
                if results.(lastKey1).rows > 1 && results.(lastKey1).cols > 1
                    results.(lastKey1).(lastKey2) = reshape(results.(lastKey1).(lastKey2),results.(lastKey1).cols,results.(lastKey1).rows)';
                end
            end
            
% %             if length(results.(lastKey1).(lastKey2)) == 9
% %                 
% %             end
        end
    end
    
    % get the value, ignoring any comments (remove whitespace)
    value = strsplit(thisLine(sepIndex+1:end), '#');
    value = strtrim(value{1}); 
    
    % attempt to convert value to numeric type
    [convertedValue, success] = str2num(value);
    if success
        value = convertedValue;
    end
    
    % store the key and value in the results
    if ~isempty(key)
%         results.(key) = value;
        if strcmpi(value(1), '!')
           setProp = 1;
           lastKey1 = key;
%            results.(key) = value;
        end
        if strcmpi(value(1),'[')
            catData = 1;
            lastKey2 = key;
            results.(lastKey1).(key) = value;
        else
            catData = 0;
            setProp = 0;
        end
        if length(key) == 4
            if strcmpi(key(1:4),'rows') || strcmpi(key(1:4),'cols')
                results.(lastKey1).(key) = round((value));
            end
        end
    end
end

end
