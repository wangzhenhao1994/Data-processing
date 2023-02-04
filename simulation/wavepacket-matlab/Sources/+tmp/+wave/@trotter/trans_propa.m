%----------------------------------------------------------------
% Perform propagation associated with transition dipoles (matrix)
%----------------------------------------------------------------
function trans_propa ( obj, psi, e, recalc)
global hamilt time
persistent dip_diagonal

if isfield(hamilt, 'dip')
    for p = 1:length(hamilt.dip)
        if ~isempty (hamilt.dip{p})
            if hamilt.dip_trans(p)
                if abs(e(p))>0
                    
                    if recalc  % || isempty(diagonal_x)
                        
                        % HACK: If we allow complex electric fields (i.e., rotating wave
                        % approximation), the eigenvectors become dependant on the phase of the
                        % electric field, and we need to recalculate them over and over again.
                        if time.efield.complex
                            
                            prt.error ('Code missing here')
                            
                            %             phase = e_x / abs(e_x);
                            %
                            %             hamilt.d_x.eig_vecs{1,1} = -ones(size(hamilt.pot{1,1}.grid)) * phase / sqrt(2);
                            %             hamilt.d_x.eig_vecs{2,1} =  ones(size(hamilt.pot{1,1}.grid))         / sqrt(2);
                            %             hamilt.d_x.eig_vecs{1,2} =  ones(size(hamilt.pot{1,1}.grid)) * phase / sqrt(2);
                            %             hamilt.d_x.eig_vecs{2,2} =  ones(size(hamilt.pot{1,1}.grid))         / sqrt(2);
                            %
                            %             for m = 1:hamilt.coupling.n_eqs
                            %                 diagonal_x{m} = exp(-abs(e_x)*hamilt.d_x.eig_vals{m});
                            %             end
                        else
                            % Profiling has shown that for more realistic settings (H2+ in a laser
                            % field with 3 coupled surfaces) the exponentiation here easily makes up
                            % 1/4-1/3 of the computing time. By calculating it only once per
                            % timestep, we can thus shave off >10%.
                            for m = 1:hamilt.coupling.n_eqs
                                dip_diagonal{p}{m} = exp( - e(p) * hamilt.dip_eig_vals{p}{m});
                            end
                        end
                    end
                    
    
                    % Transform to adiabatic (transformation matrix D+)
                    for m = 1:hamilt.coupling.n_eqs
                        psi.new{m} = zeros(size(psi.dvr{1}));
                        for n=1:hamilt.coupling.n_eqs
                            psi.new{m} = psi.new{m} ...
                                + conj(hamilt.dip_eig_vecs{p}{n,m}) .* psi.dvr{n};
                        end
                        
                        % Propagate in adiabatic representation
                        psi.new{m} = psi.new{m} .* dip_diagonal{p}{m};
                    end
    
                    % Transform back to diabatic (transformation matrix D)
                    for m = 1:hamilt.coupling.n_eqs
                        psi.dvr{m} = zeros(size(psi.new{1}));
                        for n=1:hamilt.coupling.n_eqs
                            psi.dvr{m} = psi.dvr{m} ...
                                + hamilt.dip_eig_vecs{p}{m,n} .* psi.new{n};
                        end
                    end
    
                end
            end
        end
    end
end

end
            
   