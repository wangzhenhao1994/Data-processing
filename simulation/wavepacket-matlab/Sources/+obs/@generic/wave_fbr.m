% Evaluate FBR ("momentum") expectation values and their uncertainties
function wave_fbr (obj, psi, step)
global hamilt space

if ~strcmpi( obj.qua,'mom')
    prt.error ('Invalid choice of quantity for mean value / uncertainty')
end

P = zeros(hamilt.coupling.n_eqs,1);
for m = 1:hamilt.coupling.n_eqs
    
    % If population exceeds pre-specified threshold
    P(m)  = sum ( abs(psi.dvr{m}(:)).^2 .* space.weight(:) );
    if P(m) > eps
        
        % Now integrate in the spectral basis to get the "momentum" expectation values.
        % momentum applies the momentum operator to the input.
        psi.mom = momentum(space.dof{obj.ind}, psi.dvr{m});
        Q1 = sum ( conj(psi.dvr{m}(:)) .* psi.mom(:) .* space.weight(:) ) / P(m);
        
        psi.mom = momentum(space.dof{obj.ind}, psi.mom);
        Q2 = sum ( conj(psi.dvr{m}(:)) .* psi.mom(:) .* space.weight(:) ) / P(m);
        
        % Remove the artificial(!) imaginary part here (should be tiny)
        obj.cha{m}(step) = math.real (             Q1     );
        obj.unc{m}(step) = math.real ( sqrt ( Q2 - Q1^2 ) );
        
    end
    
end
end