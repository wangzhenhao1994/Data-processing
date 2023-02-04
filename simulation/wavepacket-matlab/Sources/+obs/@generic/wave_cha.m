% Expectation values / uncertainties for each CHANNEL
function wave_cha (obj, psi, step)
global hamilt space

P = zeros(hamilt.coupling.n_eqs,1);
for m = 1:hamilt.coupling.n_eqs
    
    % If population exceeds pre-specified threshold
    P(m)  = sum ( abs(psi.dvr{m}(:)).^2 .* space.weight(:) );
    if P(m) > eps
        
        switch lower(obj.qua)
            case 'amo'
                if isempty(hamilt.amo{obj.ind}.dvr)
                    Q = 0;
                else
                    Q = hamilt.amo{obj.ind}.dvr;
                end
            case 'pos'
                Q = space.dvr{obj.ind};
            case 'pot'
                Q = hamilt.pot{m,m}.dvr;
            otherwise
                prt.error ('Invalid choice of quantity for mean value / uncertainty')
        end
        if ~isempty(Q)
            Q1 = sum (  Q(:)    .*abs(psi.dvr{m}(:)).^2 .* space.weight(:) ) / P(m);
            Q2 = sum (  Q(:).^2 .*abs(psi.dvr{m}(:)).^2 .* space.weight(:) ) / P(m);
        else
            Q1 = 0;
            Q2 = 0;
        end
        
        % Remove the artificial(!) imaginary part here (should be tiny)
        obj.cha{m}(step) = math.real (             Q1     );
        obj.unc{m}(step) = math.real ( sqrt ( Q2 - Q1^2 ) );
    end
end
end
