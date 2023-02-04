% Dichotomy gradient descent method
% probable is of same dimension as ind_c
function probable = dich_gd (obj , ind_m , ind_c , m , n )
global space time
       
    pos_new    = cell(space.n_dim,1);
    mom_new    = cell(space.n_dim,1);
    pos_old    = cell(space.n_dim,1);
    mom_old    = cell(space.n_dim,1);
    mom_oldold = cell(space.n_dim,1);
    
    frc_m_new       = cell(space.n_dim,1);
    frc_m_old       = cell(space.n_dim,1);
    frc_mats        = cell(space.n_dim,1);
    frc_mats_old    = cell(space.n_dim,1);
    frc_mats_oldold = cell(space.n_dim,1);

    pot_mats        = obj.pot_mat_old(:,:,ind_m(ind_c));
    pot_mats_old    = obj.pot_mat_oldold(:,:,ind_m(ind_c));
    pot_mats_oldold = obj.pot_mat(:,:,ind_m(ind_c));
    
    for d=1:space.n_dim
        pos_new{d} = obj.pos_old{d}(ind_m(ind_c));
        mom_new{d} = obj.mom_old{d}(ind_m(ind_c));
        
        pos_old{d} = obj.pos_oldold{d}(ind_m(ind_c));
        mom_old{d} = obj.mom_oldold{d}(ind_m(ind_c));
        
        mom_oldold{d} = obj.mom{d}(ind_m(ind_c));
        
        frc_m_new{d}  = obj.frc_old{d}(ind_m(ind_c));
        frc_m_old{d}  = obj.frc_oldold{d}(ind_m(ind_c));
        
        frc_mats{d}        = obj.frc_mat_old{d}(:,:,ind_m(ind_c));
        frc_mats_old{d}    = obj.frc_mat_oldold{d}(:,:,ind_m(ind_c));
        frc_mats_oldold{d} = obj.frc_mat{d}(:,:,ind_m(ind_c));
    end

    U_new    = obj.U_old(:,:,ind_m(ind_c));
    U_old    = obj.U_oldold(:,:,ind_m(ind_c));
    U_oldold = obj.U_new(:,:,ind_m(ind_c));
    
    D_new    = obj.D_old(:,ind_m(ind_c));
    D_old    = obj.D_oldold(:,ind_m(ind_c));
    D_oldold = obj.D_new(:,ind_m(ind_c));
    
    % back up for momentum rescaling
    frc_n_backup = ham.frc_adi ( pot_mats, frc_mats , U_new , n );
    frc_n_new    = frc_n_backup;
    
    % First derivative of gap. Computing with exact formula.
    gap_new_d1 = zeros(length(ind_c),1);
    gap_old_d1 = zeros(length(ind_c),1);
    for d = 1:space.n_dim
        gap_new_d1(:) = gap_new_d1(:) + (-1)^(n<m) * (frc_m_new{d} - frc_n_new{d}) .* mom_new{d} / space.dof{d}.mass;
    end
    
    % current gap and its second derivative
    gap_new    = abs(D_new(m,:)   -D_new(n,:)   )';
    gap_old    = abs(D_old(m,:)   -D_old(n,:)   )';
    gap_oldold = abs(D_oldold(m,:)-D_oldold(n,:))';
    
    % sum of all t_delta used for the gradient method
    sum_delta_t = zeros(length(ind_c),1);
    
    % current time step
    t_delta     = - sign(gap_new_d1) .* time.steps.s_delta * 0.5;
    t_delta_old = time.steps.s_delta * ones(length(ind_c),1);
    
    % For LZGD variant of lz_2:
    % Local minimum adiabatic energy gap (previous step) of the trajectories
    % with second time derivatives
    % https://en.wikipedia.org/wiki/Finite_difference_coefficient
    % Quite confusing since gap_oldold is the current gap!
    gap_new_d2  = (gap_oldold - 2* gap_new + gap_old) / time.steps.s_delta^2;
    
    prob_new    = lz_formula(obj, mom_new,  pot_mats, frc_mats, ...
                            U_new, D_new, gap_new, gap_new_d2, m, n); % current lz probability
    prob_old    = lz_formula(obj, mom_old,  pot_mats_old, frc_mats_old, ...
                            U_old, D_old, gap_old, gap_new_d2, m, n); % init for the gradient descent method: use same gap_new_d2
    prob_oldold = lz_formula(obj, mom_oldold,  pot_mats_oldold, frc_mats_oldold, ...
                            U_oldold, D_oldold, gap_oldold, gap_new_d2, m, n); % init for the gradient descent method: use same gap_new_d2

    for k = 1:obj.max_round_gd
        
        % Convergence criterion for gradient descent
        acc = obj.acc_gd;
        ind_k = ~ (   abs(prob_new - prob_old)    < acc * prob_old    ...
                    & abs(prob_new - prob_oldold) < acc * prob_oldold ...
                  );
        n_k = sum(ind_k);
        if (n_k == 0)
            break
        end
        
        % Interval halving method or dichotomy method
        t_delta(ind_k)    = - sign(gap_new_d1(ind_k)) .* time.steps.s_delta * 0.5^k;
        
        % sum of all t_delta used for the gradient method
        sum_delta_t(ind_k) = sum_delta_t(ind_k) + t_delta(ind_k);
        
        % position integrator
        % J. Dummer modified:
        pos_current = pos_new;
        for d = 1:space.n_dim
            pos_new{d}(ind_k) = pos_new{d}(ind_k) + (pos_new{d}(ind_k)-pos_old{d}(ind_k)) .* t_delta(ind_k) ./ t_delta_old(ind_k) + ...
                            t_delta(ind_k) .* frc_m_new{d}(ind_k) / space.dof{d}.mass .* (t_delta(ind_k) + t_delta_old(ind_k))/2;
        end
        pos_old = pos_current;
        
        pos_new_ind_k = cell(space.n_dim,1);
        for d = 1:space.n_dim
            pos_new_ind_k{d} = pos_new{d}(ind_k);
        end
        
        % Calculate and save all diabatic potentials in advance
        pot_mats_ind_k = ham.pot_dia(pos_new_ind_k);
        pot_mats(:,:,ind_k)  = pot_mats_ind_k;
        
        % Compute adiabatic potential matrix and eigenvector matrix
        [U_new(:,:,ind_k) , D_new(:,ind_k)] = ham.pot_eig_adi(pot_mats_ind_k);
        
        % Calculate and save all diabatic forces in advance
        frc_mats_ind_k = ham.frc_dia(pos_new_ind_k);
        
        frc_m_oldold = frc_m_old;
        frc_m_old   = frc_m_new;
        
        frc_m_new_ind_k = ham.frc_adi ( pot_mats_ind_k, frc_mats_ind_k , U_new(:,:,ind_k) , m );
        frc_n_new_ind_k = ham.frc_adi ( pot_mats_ind_k, frc_mats_ind_k , U_new(:,:,ind_k) , n );
        
        for d = 1:space.n_dim
            frc_m_new{d}(ind_k) = frc_m_new_ind_k{d};
            frc_n_new{d}(ind_k) = frc_n_new_ind_k{d};
            
            frc_mats{d}(:,:,ind_k) = frc_mats_ind_k{d};
        end
        
        % momentum integrator
        % J. Dummer modified:
        mom_current = mom_new;
        mom_new_ind_k = cell(space.n_dim,1);
        for d = 1:space.n_dim
            mom_new{d}(ind_k) = mom_new{d}(ind_k) + (mom_new{d}(ind_k)-mom_old{d}(ind_k)) .* t_delta(ind_k) ./ t_delta_old(ind_k) ...
                                    + t_delta(ind_k) .* (frc_m_new{d}(ind_k)-frc_m_oldold{d}(ind_k))/2;
            mom_new_ind_k{d} = mom_new{d}(ind_k);
        end
        mom_old = mom_current;

        % new gap
        gap_new(ind_k) = abs( D_new(n,ind_k) - D_new(m,ind_k) )';

        % new first derivative of gap
        % exact formula
        gap_old_d1(ind_k) = gap_new_d1(ind_k);
        gap_new_d1(ind_k) = zeros(n_k,1);
        for d = 1:space.n_dim
            gap_new_d1(ind_k) = gap_new_d1(ind_k) + (-1)^(n<m) * (frc_m_new{d}(ind_k) - frc_n_new{d}(ind_k)) .* mom_new{d}(ind_k) / space.dof{d}.mass;
        end

        % save current time step 
        t_delta_old(ind_k) = t_delta(ind_k);
        
        prob_oldold(ind_k) = prob_old(ind_k);
        prob_old(ind_k) = prob_new(ind_k);
        
        % For LZGD variant of lz_2: new second derivative of gap
        % first order finite difference
        gap_new_d2(ind_k) = (gap_new_d1(ind_k) - gap_old_d1(ind_k)) ./ t_delta(ind_k);
        
        prob_new(ind_k) = lz_formula(obj, mom_new_ind_k,  pot_mats_ind_k, frc_mats_ind_k, ...
                            U_new(:,:,ind_k), D_new(:,ind_k), gap_new(ind_k), gap_new_d2(ind_k), m, n);
        
        % Counting the extra diagonalizations produced by the gradient descent method
        time.counter.diag.gd = time.counter.diag.gd + n_k;
        
    end
    
    % Update only the probability of those trajectories which energy gap has been reduced
    probable = prob_new;
    
    
    % For momentum rescaling: Save last values of the gradient descent method.
    obj.hop_delta_t  {m,n}(ind_m(ind_c)) = time.steps.s_delta - sum_delta_t;
    obj.hop_interval{m,n}(ind_m(ind_c)) = abs(t_delta);
    
    obj.hop_pot_mat{m,n}(:,:,ind_m(ind_c)) = pot_mats;
    obj.hop_U{m,n}(:,:,ind_m(ind_c))       = U_new;
    obj.hop_D{m,n}(:,ind_m(ind_c))         = D_new;
    
    obj.hop_pos{m,n} = obj.pos_old;
    obj.hop_mom{m,n} = obj.mom_old;
    for d = 1:space.n_dim
        
        obj.hop_pos{m,n}{d}(ind_m(ind_c)) = pos_new{d};
        obj.hop_mom{m,n}{d}(ind_m(ind_c)) = mom_new{d};
        
        obj.hop_frc_n{m,n}{d}(ind_m(ind_c)) = frc_n_new{d};
        
        obj.hop_frc_mat{m,n}{d}(:,:,ind_m(ind_c)) = frc_mats{d};
    end
    
    obj.hop_ham{m,n}{m,m}(ind_m(ind_c)) = D_new(m,:);
    obj.hop_ham{m,n}{n,n}(ind_m(ind_c)) = D_new(n,:);

end

