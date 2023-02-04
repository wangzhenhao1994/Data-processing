%----------------------------------------------------------------
% Initialize propagator associated with potential energy (matrix)
%----------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function pot_init (obj)
global hamilt space time

% Detect whether any potential coupling exists
hamilt.pot_couple = false;
for m=1:hamilt.coupling.n_eqs
    for n=m+1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.pot{m,n}.dvr)
            hamilt.pot_couple = true;
        end
    end
end

% Diabatic potential coupling is absent
if  ~hamilt.pot_couple
    % prt.disp(['Initialize diagonal potential propagator for time step fraction: ', num2str(obj.frac)])
    for m=1:hamilt.coupling.n_eqs
        if ~isempty (hamilt.pot{m,m}.dvr)
            hamilt.pot_expo{m,m} = exp ( -1i * time.steps.s_delta*obj.frac * hamilt.pot{m,m}.dvr );
        else
            hamilt.pot_expo{m,m} = [];
        end
    end
    
    % Potential coupling exists
else
    
    switch hamilt.coupling.n_eqs
        
        % Scalar exponential for one equation (nothing to be done)
        case 1
            
            % Analytic matrix exponential for two equations
            % expm(-i*tau*[eta+delta beta; beta eta-delta])
        case 2
            % prt.disp(['Initialize analytical potential coupling propagator for time step fraction: ', num2str(obj.frac)])
            delta = ( hamilt.pot{1,1}.dvr - hamilt.pot{2,2}.dvr ) /2;
            eta   = ( hamilt.pot{1,1}.dvr + hamilt.pot{2,2}.dvr ) /2;
            beta  =   hamilt.pot{1,2}.dvr;
            rho   = sqrt (delta.^2 + beta.^2);
            
            ex = exp( -1i * time.steps.s_delta*obj.frac * eta );
            co = cos(       time.steps.s_delta*obj.frac * rho );
            si = sin(       time.steps.s_delta*obj.frac * rho );
            
            hamilt.pot_expo{1,1}  = ex .* (co - 1i*delta.*si./rho);
            hamilt.pot_expo{2,2}  = ex .* (co + 1i*delta.*si./rho);
            hamilt.pot_expo{1,2}  = ex .* (   - 1i* beta.*si./rho);
            
            % Special care needs to be taken at conical interaction(s)
            % alfa=gamma=eta, delta=0, beta=0, rho=0
            % Hence, potential matrix is diagonal
            ci = find(rho(:)==0); % Using linear indexing
            hamilt.pot_expo{1,1}(ci)  = ex(ci);
            hamilt.pot_expo{2,2}(ci)  = ex(ci);
            hamilt.pot_expo{1,2}(ci)  = 0;
            
            % Numerical matrix exponential for more than two coupled equations
        otherwise
            
            % prt.disp(['Initialize numerical potential coupling propagator for time step fraction: ', num2str(obj.frac)])
                        
            % Loop over all points of (multidimensional) grid
            for gr=1:space.n_tot
                
                % Preallocate (should become 3D array, for use with "reshape")
                pot_mat = zeros(hamilt.coupling.n_eqs);
                
                % Extract diabatic potentials from cell array
                for m=1:hamilt.coupling.n_eqs
                    if ~isempty(hamilt.pot{m,m}.dvr)
                        pot_mat(m,m) = hamilt.pot{m,m}.dvr(gr);
                    end
                        
                    % Diabatic potential coupling (Hermitian)
                    for n=m+1:hamilt.coupling.n_eqs
                        if ~isempty(hamilt.pot{m,n}.dvr)
                            pot_mat(m,n) = hamilt.pot{m,n}.dvr(gr);
                            pot_mat(n,m) = conj(pot_mat(m,n));
                        end
                    end
                    
                end
                
                % Either use Matlab's expm function ...
                % exp_pot = expm ( -i * time.steps.s_delta*obj.frac * pot_mat );
                
                % ... or use eigen/values/vectors (from expmdemo3.m)
                [V,D] = eig(pot_mat);
                exp_pot = V * diag(exp(-1i * time.steps.s_delta * obj.frac * diag(D))) / V;
                
                % Store potential propagator for later use
                % Real, symmetic matrices
                for m=1:hamilt.coupling.n_eqs
                    for n=m:hamilt.coupling.n_eqs
                        hamilt.pot_expo{m,n}(gr,1) = exp_pot(m,n);
                    end
                end
                
            end % loop over grid points
            
            % If necessary, reshape vectors into arrays
            if space.n_dim>1
                n_pts = zeros (1, space.n_dim); % row vector
                for dim = 1:space.n_dim
                    n_pts(dim) = space.dof{dim}.n_pts;
                end
                
                for m=1:hamilt.coupling.n_eqs
                    for n=m:hamilt.coupling.n_eqs
                        hamilt.pot_expo{m,n} = reshape ( hamilt.pot_expo{m,n}, n_pts );
                    end
                end
            end
            
    end % switch n_eqs
    
end % if potential coupling exists

end
        