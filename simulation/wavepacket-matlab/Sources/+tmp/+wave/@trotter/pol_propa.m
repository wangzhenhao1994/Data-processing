            
%---------------------------------------------------------------
% Perform propagation associated with polarizabilities (vector)
%---------------------------------------------------------------
function pol_propa ( obj,psi, e )
global hamilt

if isfield(hamilt, 'pol')
    for p = 1:size(hamilt.pol,1)
        if abs(e(p))>0
            for q = p:size(hamilt.pol,2)
                if abs(e(q))>0
                    if ~isempty (hamilt.pol{p,q})
                        for m = 1:hamilt.coupling.n_eqs
                            if  ~isempty ( hamilt.pol_prod{p,q}{m} )
                                psi.dvr{m} = exp(hamilt.pol_prod{p,q}{m} * (-e(p)*e(q)/2)) .* psi.dvr{m};
                            end
                        end
                    end
                end
            end
        end
    end 
end

end
