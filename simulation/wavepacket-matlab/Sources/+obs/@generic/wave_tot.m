% Evaluate expectation values for TOTAL wave function
function wave_tot (obj, psi, step)
global hamilt space

total = 0;
P = zeros(hamilt.coupling.n_eqs,1);
for m = 1:hamilt.coupling.n_eqs
    
    % If population exceeds pre-specified threshold
    P(m)  = sum ( abs(psi.dvr{m}(:)).^2 .* space.weight(:) );
    if P(m) > eps
        if ~strcmpi (obj.qua,'pop')
            total = total + obj.cha{m}(step) * P(m);
        else
            total = total + obj.cha{m}(step);
        end
    end
end

if ~strcmpi (obj.qua,'pop')
    total = total / sum(P);
end

% Remove the artificial(!) imaginary part here (should be tiny)
obj.tot (step) = math.real (total);

end
 