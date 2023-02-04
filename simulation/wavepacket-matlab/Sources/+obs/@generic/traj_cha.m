% Evaluate expectation values / uncertainties for each CHANNEL
function traj_cha (obj, state, step)
global hamilt

for m = 1:hamilt.coupling.n_eqs
    
    % If population exceeds pre-specified threshold
    N  = nnz ( state.cha==m );
    if N
        
        switch lower(obj.qua)
            case 'amo'
                if isempty(hamilt.amo{obj.ind}.dvr) % ?????? Scheisse
                    Q = 0;
                else
                    Q = hamilt.amo{obj.ind}.dvr; % ????? Scheisse
                end
            case 'pos'
                Q = state.pos{obj.ind}(state.cha==m);
            case 'mom'
                Q = state.mom{obj.ind}(state.cha==m);
            case 'pot'
                Q = state.pot(state.cha==m);
            case 'kin'
                Q = state.kin(state.cha==m);
            otherwise
                prt.error ('Invalid choice of quantity for mean value / uncertainty')
        end
        Q1 = sum (  Q(:)     ) / N;
        Q2 = sum (  Q(:).^2  ) / N;
        
        % Remove the artificial(!) imaginary part here (should be tiny)
        obj.cha{m}(step) = math.real (             Q1     );
        obj.unc{m}(step) = math.real ( sqrt ( Q2 - Q1^2 ) );
    end
end
end