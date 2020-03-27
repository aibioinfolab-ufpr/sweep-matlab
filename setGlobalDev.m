function setGlobalDev(val)
% The global DEV_MODE variable must be seted to test the complete 
% application as the compiled version, otherwise the default mode to run
% in the IDE is without file generation.
global DEV_MODE
DEV_MODE = val; end