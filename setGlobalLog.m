function setGlobalLog(val)
% The global GLOBAL_LOG variable is the only seted outside the main sweep 
% script, so it requires an individual set function file.
global GLOBAL_LOG
GLOBAL_LOG = val; end