% Evaluate kinetic energy expectation values and their uncertainties
function wave_kin (obj, psi, step)
global hamilt space

if ~strcmpi( obj.qua,'kin')
    prt.error ('Invalid choice of quantity for mean value / uncertainty')
end

Q1 = zeros(hamilt.coupling.n_eqs,1);
Q2 = zeros(hamilt.coupling.n_eqs,1);
P  = zeros(hamilt.coupling.n_eqs,1);
for m = 1:hamilt.coupling.n_eqs
    P(m)  = sum ( abs(psi.dvr{m}(:)).^2 .* space.weight(:) );
end

for k = 1:space.n_dim
    
    % Calculate T|Psi> and put it in psi.newn and get <E_kin>
    kinetic(space.dof{k}, psi, false);
    for m = 1:hamilt.coupling.n_eqs
        if P(m) > eps % If population exceeds pre-specified threshold
            Q1(m) = Q1(m) + sum( conj(psi.new{m}(:)) .* psi.dvr{m}(:) .* space.weight(:)) / P(m);
        end
    end
    
    % Apply T another time on psi to get <E_kin^2>
    kinetic(space.dof{k}, psi, true);
    for m = 1:hamilt.coupling.n_eqs
        if P(m) > eps % If population exceeds pre-specified threshold
            Q2(m) = Q2(m) + sum( conj(psi.new{m}(:)) .* psi.dvr{m}(:) .* space.weight(:)) / P(m);
        end
    end
end

% Now if T = T1 + T2 + ..., <E_kin2 ^> has also cross terms
for k = 1:space.n_dim
    for l = 1:space.n_dim
        if k == l
            continue
        end
        kinetic(space.dof{k}, psi, false);
        kinetic(space.dof{l}, psi, true);
        for m = 1:hamilt.coupling.n_eqs
            if P(m) > eps
                Q2(m)  = Q2(m) + ...
                    sum( conj(psi.new{m}(:)) .* psi.dvr{m}(:) ...
                    .* space.weight(:)) / P(m);
            end
        end
    end
end

% For extra kinetic energy operators, the same thing has to be done as well.
if isfield(hamilt, 'kin')
    for n = 1:length(hamilt.kin)
        
        % Calculate T|Psi> and put it in psi.newn and get <E_kin>
        kinetic(hamilt.kin{n}, psi, false);
        for m = 1:hamilt.coupling.n_eqs
            if P(m) > eps
                Q1(m) = Q1(m) + sum( conj(psi.new{m}(:)) .* psi.dvr{m}(:) .* space.weight(:)) / P(m);
            end
        end
        
        kinetic(hamilt.kin{n}, psi, true);
        for m = 1:hamilt.coupling.n_eqs
            if P(m) > eps
                Q2(m) = Q2(m) + sum( conj(psi.new{m}(:)) .* psi.dvr{m}(:) .* space.weight(:)) / P(m);
            end
        end
    end
    
    % Cross-terms. (1) grid+external, external+grid,
    for n = 1:length(hamilt.kin)
        for k = 1:space.n_dim
            kinetic(hamilt.kin{n},psi, false);
            kinetic(space.dof{k}, psi, true);
            for m = 1:hamilt.coupling.n_eqs
                if P(m) > eps
                    Q2(m)  = Q2(m) + sum( conj(psi.new{m}(:)) .* psi.dvr{m}(:) .* space.weight(:)) / P(m);
                end
            end
            
            % Cross-terms. (2) external+grid,
            kinetic(space.dof{k}, psi, false);
            kinetic(hamilt.kin{n},psi, true);
            for m = 1:hamilt.coupling.n_eqs
                if P(m) > eps
                    Q2(m)  = Q2(m) + sum( conj(psi.new{m}(:)) .* psi.dvr{m}(:) .* space.weight(:)) / P(m);
                end
            end
        end
    end
    
    % Cross-terms. (3) external+external
    for n = 1:length(hamilt.kin)
        for o = 1:length(hamilt.kin)
            if n == o
                continue;
            end
            kinetic(hamilt.kin{n}, psi, false);
            kinetic(hamilt.kin{n}, psi, true);
            for m = 1:hamilt.coupling.n_eqs
                if P(m) > eps
                    Q2(m)  = Q2(m) + sum( conj(psi.new{m}(:)) .* psi.dvr{m}(:) .* space.weight(:)) / P(m);
                end
            end
        end
    end
end

% Finally, save the result
for m = 1:hamilt.coupling.n_eqs
    
    % If population exceeds pre-specified threshold
    if P(m) > eps
        
        % Remove the artificial(!) imaginary part here (should be tiny)
        obj.cha{m}(step) = math.real (                Q1(m)     );
        obj.unc{m}(step) = math.real ( sqrt ( Q2(m) - Q1(m)^2 ) );
    end
end
end
