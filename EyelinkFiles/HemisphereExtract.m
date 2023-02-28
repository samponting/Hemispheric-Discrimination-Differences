% Define the filename prefixes to search for
filenamePrefixes = {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7'};

% Initialize an empty structure to store the matched strings for each prefix
HemispheresByPPT = struct();

% Loop through each filename prefix
for i = 1:length(filenamePrefixes)
    
    % Get a list of files that match the prefix
    prefix = filenamePrefixes{i};
    files = dir([prefix, '*']);
    
    % Initialize an empty structure to store the matched strings for this prefix
    prefixMatchedStrings = struct();
    
    % Loop through each file
    for j = 1:length(files)
        
        % Load the mat file into a variable called matFile
        load(files(j).name, 'matFile');
        
        % Initialize an empty cell array to store the matched strings
        matchedStrings = {};
        
        % Loop through each cell in the info array
        for k = 1:length(matFile.Events.Messages.info)

            % Check if the cell contains the string 'Hemisphere: Lower' or 'Hemisphere: Upper'
            if contains(matFile.Events.Messages.info{k}, 'Hemisphere: Lower') || contains(matFile.Events.Messages.info{k}, 'Hemisphere: Upper')

                % Add the matched string to the end of the matchedStrings array
                matchedStrings{end+1} = matFile.Events.Messages.info{k};

            end

        end
        
        % Store the matched strings in the prefixMatchedStrings structure using the filename as the fieldname
        [~, name, ~] = fileparts(files(j).name);
        prefixMatchedStrings.(name) = matchedStrings;
        
    end
    
    % Store the prefixMatchedStrings structure in the matchedStringsByPrefix structure using the prefix as the fieldname
    HemispheresByPPT.(prefix) = prefixMatchedStrings;
    
end

% Display the matched strings for each prefix in the command window
disp(HemispheresByPPT);

save('Hemisphere_Extract.mat','HemispheresByPPT')
