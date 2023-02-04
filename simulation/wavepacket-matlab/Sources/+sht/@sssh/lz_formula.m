%--------------------------------------------------------------------------
%
% Landau-Zener (LZ) based transition probabilities for use
% in "single switch surface hopping" trajectory simulations
%
% Several different LZ variants are available
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function probable = lz_formula(obj, mom, dia_pot_mats, dia_frc_mats, U, D, gap, gap_d2, m, n)
global hamilt space

switch obj.lz_variant
    case 1
        % Eq. (3) from doi:10.1103/PhysRevA.84.014701 or 5 lines above
        % Eq. (3) from doi:10.1063/1.4882073 but why is gap squared there?!?
        
        nac_mn = ham.nac_mn (dia_pot_mats,dia_frc_mats,U,D,m,n);
        
        % Non-adiabatic coupling
        coupling = zeros (length(mom{1}),1);
        for d = 1:space.n_dim
            coupling = coupling + nac_mn{d} .* mom{d} / space.dof{d}.mass;
        end
        
        % Adiabatic Landau Zener probability: only for "critical" trajectories
        probable = exp ( - pi/4 * gap ./ abs(coupling) );
        
    case 2
        % Eq. (2) from doi:10.1103/PhysRevA.84.014701
        % Eq. (4) from doi:10.1063/1.4882073
        
        % Local minimum adiabatic energy gap (previous step) of the trajectories
        % with second time derivatives
        % Adiabatic Landau Zener probability: only for "critical" trajectories
        probable = exp ( - pi/2 * sqrt(gap.^3 ./ gap_d2 ) );
        
    case 3
        % Eq. (1) from doi:10.1103/PhysRevA.84.014701
        
        gap_dia_d1 = zeros (length(mom{1}),1);
        for d = 1:space.n_dim
            gap_frc    = zeros (length(mom{1}),1);
            gap_frc(:) = dia_frc_mats{d}(m,m,:) - dia_frc_mats{d}(n,n,:);
            gap_dia_d1 = gap_dia_d1 + gap_frc .* mom{d} / space.dof{d}.mass;
        end
        
        coupling    = zeros (length(mom{1}),1);
        coupling(:) = dia_pot_mats(m,n,:);
        
        probable = exp ( - 2*pi * coupling.^2 ./ abs(gap_dia_d1) );
        
    case 4
        % Eq. (3) from doi:10.1063/1.4882073
        
        coupling_1 = zeros (length(mom{1}),1);
        coupling_2 = zeros (length(mom{1}),1);
        
        for d = 1:space.n_dim
            gap_frc    = zeros (length(mom{1}),1);
            gap_frc(:) = dia_frc_mats{d}(m,m,:) - dia_frc_mats{d}(n,n,:);
            coupling_1 = coupling_1 + 0.5 * gap_frc .* mom{d} / space.dof{d}.mass;
            
            frc_mn_d    = zeros (length(mom{1}),1);
            frc_mn_d(:) = dia_frc_mats{d}(m,n,:);
            coupling_2  = coupling_2 + frc_mn_d .* mom{d} / space.dof{d}.mass;
        end
        
        coupling = sqrt(coupling_1.^2 + coupling_2.^2);
        
        probable = exp ( - pi/4 * gap.^2 ./ coupling );
        
    case {5, 6}
        % See ticket #224 1. and 2. formula
        % Currently restricted only to one single example
        
        epsilon = hamilt.pot{1,1}.eps;
        %                     epsilon = 1 / sqrt(space.dof{1}.mass);
        
        delta   = zeros (length(mom{1}),1);
        p       = zeros (length(mom{1}),1);
        
        %                     delta(:) = dia_pot_mats(m,n,:) .* epsilon;
        delta(:) = hamilt.pot{1,1}.delta;
        
        p(:) = mom{1} .* epsilon;
        %                     p(:) = mom{1};
        
        n = sign(p) .* sqrt( p.^2 - 4 * delta );
        
        t = pi * ( sqrt( 1 + delta.^2 ) - 1 );
        
        probable = exp ( - t .* abs( p - n ) ./ ( delta .* epsilon ) );
        
        if(obj.lz_variant == 6)
            
            C = ( (n+p) ./ (2 * abs(n)) ).^2;
            
            probable = C .* probable;
            
        end
        
end