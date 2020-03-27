function mret = varginfind(vin, xname, xelse)
% Retorna valor atribuido a xname se este ocorrer em vin (varargin externa)
mret = xelse;
%
if ~isempty(vin)
    for ii=1:length(vin)
        V = vin{ii};
        if ischar(V) 
            if strcmp(upper(V),upper(xname))
                mret = vin{ii+1};
            end
        end
    end
end

