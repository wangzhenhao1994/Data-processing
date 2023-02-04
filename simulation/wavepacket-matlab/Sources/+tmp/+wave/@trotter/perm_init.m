%------------------------------------------------------------------
% Initialize propagation associated with permanent dipoles (vector)
%------------------------------------------------------------------
function perm_init (obj)
global hamilt time

if isfield(hamilt, 'dip')
    for p = 1:length(hamilt.dip)
        if ~isempty (hamilt.dip{p})
            for m=1:hamilt.coupling.n_eqs
                if ~isempty ( hamilt.dip{p}{m,m}.dvr )
                    hamilt.dip_prod{p}{m} = -1i * time.steps.s_delta * obj.frac * hamilt.dip{p}{m,m}.dvr;
                else
                    hamilt.dip_prod{p}{m} = [];
                end
            end
        end
    end
end

end
