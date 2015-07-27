classdef gatingmodel
    %GATINGMODEL Class for a single channel gating mechanism
    
    properties
        numstates
        numrates
        Q
        openstates
        shutstates
        idxall
        idxvary
        idxtheta
        gamma
        xi
        U
        R1
        R2
    end
    
    methods
        function mdl = gatingmodel(varargin)
            p = inputParser;
            addOptional(p,'filename','',@(x) ischar(x) && exist(x,'file'));
            parse(p,varargin{:});
            filename = p.Results.filename;
            
            if ~isempty(filename)
                try
                    [mdl.Q, mdl.openstates, mdl.shutstates, mdl.idxall, ...
                        mdl.idxvary, tmpgamma, mdl.xi, mdl.numstates, ...
                        mdl.numrates] = qmfread(filename);
                    [mdl.U,mdl.R1,mdl.R2,mdl.idxtheta,mdl.gamma] = ...
                        constrainqr(tmpgamma,mdl.idxall,mdl.idxvary);
                catch err
                end
            end
        end
        
        function mdl = load(mdl, filename)
            if ~ischar(filename) || ~exist(filename,'file')
                return
            end
            try
                [mdl.Q, mdl.openstates, mdl.shutstates, mdl.idxall, ...
                    mdl.idxvary, tmpgamma, mdl.xi, mdl.numstates, ...
                    mdl.numrates] = qmfread(filename);
                [mdl.U,mdl.R1,mdl.R2,mdl.idxtheta,mdl.gamma] = ...
                    constrainqr(tmpgamma,mdl.idxall,mdl.idxvary);
            catch err
                err
            end
        end
    end
    
end

