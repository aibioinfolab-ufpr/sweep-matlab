function generate_log(message, type)
% Generates the log for all the execution steps. Centralizes the
% information in a single component in memory, so it can be uptdate to the
% UI and dumps to disk log file.

% ARGS:
%       message: Message to update in the log;
%       type: The tipe of message - Update, Warning or Error

% RETURNS:
%       GLOBAL_LOG: Updates the global log already in memory.

% OUTPUTS:
%       Dump the message to the global log, available in the output folder.

% Mariane Goncalves Kulik (mgkulik) - 2018-nov-18
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

% if nargin > 2
%     DEV_MODE = varargin{1};
% else
%     DEV_MODE = 0;
% end
global DEV_MODE;

if ~DEV_MODE
    global USER_PATH;
    global GLOBAL_LOG;
end

if type==0
    typeTxt = '(Update)';
elseif type==1
    typeTxt = '(Warning)';
elseif type==2
    typeTxt = '(Error)';
end


if type~=3
    strMessage = cellstr(strcat(datestr(now), {' '}, typeTxt, {' '}, message));
else
    % Added a spacer option to return to the UI. It'll not be saved to the
    % log file
    spacer = repmat(char(message),1,40);
    strMessage = cellstr(spacer);
end

% Adding message to the global memory log
if ~DEV_MODE
    setGlobalLog([GLOBAL_LOG; strMessage]);
    if type~=2
        disp(char(strMessage));
    end

    if type~=3
        % Saving the message to disk. This will append the message to the exist log
        % in the folder output or create one if not found.
        fid = fopen(fullfile(USER_PATH, 'sweeplog.txt'), 'a');
        if fid == -1
            error('Cannot open log file. Please check the permissions for the "output" folder.');
        end
        fprintf(fid, '%s\r\n', char(strMessage));
        fclose(fid);
    end
end

