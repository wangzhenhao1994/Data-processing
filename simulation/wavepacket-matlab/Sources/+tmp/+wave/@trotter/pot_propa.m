%--------------------------------------------------------------
% Perform propagation associated with potential energy (matrix)
%--------------------------------------------------------------
function pot_propa (obj, psi)
global hamilt

% Dynamics along diabatic potential energy surfaces
for m = 1:hamilt.coupling.n_eqs
    if ~isempty (hamilt.pot_expo{m,m})
        psi.new{m} = hamilt.pot_expo{m,m} .* psi.dvr{m};
    else
        psi.new{m} =                         psi.dvr{m};
    end
        
end

% Potential coupling
if hamilt.pot_couple
    for m = 1:hamilt.coupling.n_eqs
        for n = m+1:hamilt.coupling.n_eqs
            if ~isempty ( hamilt.pot_expo{m,n} )
                psi.new{m} = psi.new{m} + hamilt.pot_expo{m,n} .* psi.dvr{n};
                psi.new{n} = psi.new{n} + hamilt.pot_expo{m,n} .* psi.dvr{m};
            end
        end
    end
end

% Save
psi.dvr = psi.new;
end