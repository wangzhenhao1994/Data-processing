%------------------------------------------------------------------
% Initialize propagation associated with polarizabilities (vector)
%------------------------------------------------------------------
function pol_init (obj)
global hamilt time

if isfield(hamilt, 'pol')
    for p = 1:size(hamilt.pol,1)
        for q = p:size(hamilt.pol,2)
            if ~isempty (hamilt.pol{p,q})
                for m=1:hamilt.coupling.n_eqs
                    if ~isempty ( hamilt.pol{p,q}{m,m}.dvr )
                        hamilt.pol_prod{p,q}{m} = -1i * time.steps.s_delta * obj.frac * hamilt.pol{p,q}{m,m}.dvr;
                    else
                        hamilt.pol_prod{p,q}{m} = [];
                    end
                end
            end
        end
    end
end

end