%-------------------------------------------------------------------
% Initialize propagation associated with transition dipoles (matrix)
%-------------------------------------------------------------------
function trans_init (obj)
global hamilt space time

% Detect any transition dipole moment along x|y
if isfield(hamilt, 'dip')
    for p = 1:length(hamilt.dip)
        if ~isempty (hamilt.dip{p})
            hamilt.dip_trans(p) = false;
            for m=1:hamilt.coupling.n_eqs
                for n=m+1:hamilt.coupling.n_eqs
                    if  ~isempty(hamilt.dip{p}{m,n}.dvr)
                        hamilt.dip_trans(p) = true;
                    end
                end
            end
        end
    end
end


switch hamilt.coupling.n_eqs
    
    % No transition dipoles for one channel (nothing to be done)
    case 1
        
        % Analytical eigen/values/vectors for two coupled equations
        % [V,D]=eig(-i*tau*[0 d12;d12 0])
    case 2
        
        if isfield(hamilt, 'dip')
            for p = 1:length(hamilt.dip)
                if ~isempty (hamilt.dip{p})
                    if hamilt.dip_trans(p)
                        hamilt.dip_eig_vals{p}{1} = + 1i * time.steps.s_delta*obj.frac * hamilt.dip{p}{1,2}.dvr;
                        hamilt.dip_eig_vals{p}{2} = - 1i * time.steps.s_delta*obj.frac * hamilt.dip{p}{1,2}.dvr;
                        hamilt.dip_eig_vecs{p}{1,1} = - ones(size(hamilt.pot{1,1}.dvr))/sqrt(2);
                        hamilt.dip_eig_vecs{p}{1,2} = + ones(size(hamilt.pot{1,1}.dvr))/sqrt(2);
                        hamilt.dip_eig_vecs{p}{2,1} = + ones(size(hamilt.pot{1,1}.dvr))/sqrt(2);
                        hamilt.dip_eig_vecs{p}{2,2} = + ones(size(hamilt.pot{1,1}.dvr))/sqrt(2);
                    end
                end
            end 
        end
        
        % Numerical eigen/values/vectors for more than two coupled equations
    otherwise
        
        % Transition dipole moments along x|y direction
        if isfield(hamilt, 'dip')
            for p = 1:length(hamilt.dip)
                if ~isempty (hamilt.dip{p})
                    if hamilt.dip_trans(p)
                        
                        % Preallocate
                        dip_mat{p} = zeros(hamilt.coupling.n_eqs);
                        hamilt.dip_eig_vals{p} = cell(hamilt.coupling.n_eqs,1);
                        hamilt.dip_eig_vecs{p} = cell(hamilt.coupling.n_eqs  );
                        for m=1:hamilt.coupling.n_eqs
                            hamilt.dip_eig_vals{p}{m} = zeros(size(hamilt.pot{1,1}.dvr));
                            for n=1:hamilt.coupling.n_eqs % Full matrices: not symmetric
                                hamilt.dip_eig_vecs{p}{m,n} = zeros(size(hamilt.pot{1,1}.dvr));
                            end
                        end
                                               
                        % Loop over all grid points
                        for gr=1:space.n_tot
                            
                            % Assemble matrix of available transition dipoles
                            for m=1:hamilt.coupling.n_eqs
                                for n=m+1:hamilt.coupling.n_eqs
                                    if  ~isempty(hamilt.dip{p}{m,n}.dvr)
                                        dip_mat{p}(m,n) = hamilt.dip{p}{m,n}.dvr(gr);
                                        dip_mat{p}(n,m) = hamilt.dip{p}{m,n}.dvr(gr);
                                    end
                                end % for n
                            end % for m
                            
                            % Numerical eigen/values/vectors, save for later use
                            [V,D] = eig ( dip_mat{p} );
                            for m=1:hamilt.coupling.n_eqs
                                % Don't put the factors in the eigenvector routine! It makes
                                % the calculation slower (complex instead of real arithmetic),
                                % and can lead to problems with non-orthonormalized eigenvectors!
                                hamilt.dip_eig_vals{p}{m}(gr) = D(m,m) * -1i * time.steps.s_delta * obj.frac;
                                for n=1:hamilt.coupling.n_eqs % Full matrix
                                    hamilt.dip_eig_vecs{p}{m,n}(gr) = V(m,n);
                                end % for n
                            end % for m
                            
                        end % for gr
                        
                        % Reshaping vectors into arrays
                        if space.n_dim>1
                            n_pts = zeros (1, space.n_dim); % row vector
                            for dim = 1:space.n_dim
                                n_pts(dim) = space.dof{dim}.n_pts;
                            end
                            
                            for m=1:hamilt.coupling.n_eqs
                                hamilt.dip_eig_vals{p}{m} = reshape ( hamilt.dip_eig_vals{p}{m}, n_pts );
                                for n=1:hamilt.coupling.n_eqs
                                    hamilt.dip_eig_vecs{p}{m,n} = reshape ( hamilt.dip_eig_vecs{p}{m,n}, n_pts );
                                end
                            end
                        end
                        
                    end % if hamilt.dip_trans(p)
                end
            end
        end
end % switch n_eqs

end
