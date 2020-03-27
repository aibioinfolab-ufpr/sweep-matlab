function mret = vararginnamed(vin, xname, xelse)
% Takes each of the named parameters received for the varargin variable and
% sets the value of the next informed parameter. It handles functions with
% multiple varagin entries.
%   PARAMETERS:
%       vin: Varargin entry
%       xname: Parameter name
%       xelse: Default value that should be used when the named parameter
%       isn't informed.
%   OUTPUT:
%       mret: 
% Roberto Tadeu Raittz (rraittz) - 2018-nov-18
% UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

mret = xelse;
%
if ~isempty(vin)
    for ii=1:length(vin)
        V = vin{ii};
        if ischar(V) 
            if strcmpi(V,xname)
                mret = vin{ii+1};
            end
        end
    end
end

