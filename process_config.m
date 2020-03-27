function conf = process_config(data_type)
% This function processes the config.ini information and set it's values to
% the global config variables, if available, otherwise set the default
% values checked in 16GB memory computer.

% As we access the fasta folder to read the files, we get the original
% source path, so we'll be able to return when required.
global DEV_MODE;

setMainPath(pwd);
MAIN_PATH = getMainPath;

conf = load_config;
% Loads the user path from the config file, if available
if isfield(conf, 'USER_PATH')
    setGlobalPath(conf.USER_PATH{:});
else
    setGlobalPath(fullfile(pwd, 'output'));
end

% Loads the user path from the config file, if available
if isfield(conf, 'PROJ_MAT_DEFAULT')
    setProjPath(conf.PROJ_MAT_DEFAULT{:});
else
    if strcmp(data_type, 'AA')
        setProjPath(fullfile(MAIN_PATH, 'orthMat_default212_AA_600.sweepproj'));
    else
        setProjPath(fullfile(MAIN_PATH, 'orthMat_default474_NT_600.sweepproj'));
    end
end

% Check if each other parameter to use only for servers has changed.
% Otherwise sets the values tested with memory of 16GB 

% GROWTH_RATE is how bigger the in memory fasta struct will be, so it'll be
% able to project its final size to manage the memory
if isfield(conf, 'GROWTH_RATE')
    setGrowthRate(str2double(conf.GROWTH_RATE{:}));
else
    setGrowthRate(1.63);
end

% MAX_MF_SIZE is the maximum size allowed in memory, considering 16GB. For
% bigger servers this value can be increased to avoid disk access.
if isfield(conf, 'MAX_MF_SIZE')
    setMaxMultiFastaSize(str2double(conf.MAX_MF_SIZE{:}));
else
    setMaxMultiFastaSize(1.2885e+09);
end

% MAX_SEQS_AMMOUNT is the maximum ammount of sequences in memory during the
% fas2sweep transformation. This measure has also taken in a 16GB memory
% machine and prevents swappping. May be incresed in computers with larger
% memory to process it faster.
if isfield(conf, 'MAX_SEQS')
    setMaxSeqsInMemory(str2double(conf.MAX_SEQS{:}));
else
    setMaxSeqsInMemory(2E+3);
end

USER_PATH = getGlobalPath;
%GROWTH_RATE = getGrowthRate;
%MAX_MF_SIZE = getMaxMultiFastaSize;
%MAX_SEQS = getMaxSeqsInMemory;
if ~DEV_MODE
    if ~exist(USER_PATH, 'dir')
       try
           mkdir(USER_PATH);
       catch ME
           error('Unable to create the user path where the app outputs will be dumped. Make sure the app has access to its souce folder or inform a new path in the configuration file.')
       end
    end
end

end


% Setting global constants
function setGlobalPath(val)
global USER_PATH
USER_PATH = val; end

function setGrowthRate(val)
global GROWTH_RATE
GROWTH_RATE = val; end

function setMaxMultiFastaSize(val)
global MAX_MF_SIZE
MAX_MF_SIZE = val; end

function setMaxSeqsInMemory(val)
global MAX_SEQS
MAX_SEQS = val; end

function setMainPath(val)
global MAIN_PATH
MAIN_PATH = val; end

function setProjPath(val)
global PROJ_PATH
PROJ_PATH = val; end

% Getting global constants
function x = getGlobalPath; global USER_PATH; x = USER_PATH; end
function x = getGrowthRate; global GROWTH_RATE; x = GROWTH_RATE; end
function x = getMaxMultiFastaSize; global MAX_MF_SIZE; x = MAX_MF_SIZE; end
function x = getMaxSeqsInMemory; global MAX_SEQS; x = MAX_SEQS; end
function x = getMainPath; global MAIN_PATH; x = MAIN_PATH; end
function x = getProjPath; global PROJ_PATH; x = PROJ_PATH; end