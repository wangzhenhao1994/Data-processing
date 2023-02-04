% Evaluate expectation values for TOTAL wave function
function traj_tot (obj, state, step)
global hamilt

total = 0;
for m = 1:hamilt.coupling.n_eqs
    
    % If population exceeds pre-specified threshold
    N  = nnz ( state.cha==m );
    if N
        if ~strcmpi (obj.qua,'pop')
            total = total + obj.cha{m}(step) * N;
        else
            total = total + obj.cha{m}(step);
        end
    end
end

if ~strcmpi (obj.qua,'pop')
    total = total / state.n_p;
end

% Remove the artificial(!) imaginary part here (should be tiny)
obj.tot (step) = math.real (total);

end
 