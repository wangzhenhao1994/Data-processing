%---------------------------------------------------------------
% Perform propagation associated with permanent dipoles (vector)
%---------------------------------------------------------------
function perm_propa ( obj, psi, e )
global hamilt

if isfield(hamilt, 'dip')
    for p = 1:length(hamilt.dip)
        if abs(e(p))>0
            if ~isempty (hamilt.dip{p})
                for m = 1:hamilt.coupling.n_eqs
                    if  ~isempty ( hamilt.dip_prod{p}{m} )
                        psi.dvr{m} = exp(hamilt.dip_prod{p}{m} * (-e(p))) .* psi.dvr{m};
                    end
                end
            end
        end
    end
end

end