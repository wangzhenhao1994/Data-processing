%-------------------------------------------------------------------------
%
% Transform (matrix-valued) potential and (vector-valued) wavefunction
% from diabatic to adiabatic representation or back. It is important 
% that the dia2adi transformation has to be called BEFORE the adi2dia
% because the latter one uses (diabatic) data stored by the former one.
%
% Note that within the context of fully quantum-mechanical WavePacket 
% simulations the adiabatic picture is used only for the purpose
% of visualizing thee wavefunctions and calculating corresponding 
% expectation values but NOT for the propagation itself! This is
% because of the problems associated with the (near) singularities 
% of the kinetic (derivative) couplings at (avoided) crossings and 
% intersections.
%
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008,2010 Ulf Lorenz
%
% see the README file for license details.

function adiabatic ( psi, step, direction )
global hamilt space time

% Use this method only *during* propagation
if step < 0
    return
end

% Use this method only if adiabatic representation is desired
if strcmpi(hamilt.coupling.represent,'dia')
    return
end

%% Diabatic to adiabatic ("forward transformation")
if strcmpi(direction,'dia2adi')
    
    % Instantaneous electric field present?
    if isfield(time, 'pulse') && ~isempty(time.efield.grid) 
        e = [ ...
            time.efield.grid{1}(time.steps.offset)
            time.efield.grid{2}(time.steps.offset) ];
    else
        e = [0 0];
    end
    
    % Solve eigenproblem: For time-dependent Hamiltonians
    % for every time step; else only for the first step
    if step==1 || isfield(time, 'pulse')
        
        switch hamilt.coupling.n_eqs
            case 1 % Light-dressed potential
                ham.pot_eig_1 (e, 1)
                
            case 2 % Analytic diagonalization for two coupled channels
                ham.pot_eig_2 (e, 1)
                
            otherwise % Numerical diagonalization for more than two coupled channels
                ham.pot_eig_N (e, 1)
        end
        
    end
    
    % Save diabatic potentials
    for m = 1:hamilt.coupling.n_eqs
        for n = m:hamilt.coupling.n_eqs
            if ~isempty (hamilt.pot{m,n}.dvr)
                hamilt.pot{m,n}.dia = hamilt.pot{m,n}.dvr;
            else
                hamilt.pot{m,n}.dia = [];
            end
        end
    end
    
    % Obtain adiabatic potentials: diagonal only, off-diagonals empty
    for m = 1:hamilt.coupling.n_eqs
        hamilt.pot{m,m}.dvr = hamilt.eig_val{m}; 
        for n = m+1:hamilt.coupling.n_eqs
            hamilt.pot{m,n}.dvr = [];
        end
    end
    
    % Save diabatic wave functions 
    for m = 1:hamilt.coupling.n_eqs
        psi.dia{m} = psi.dvr{m};
        psi.dvr{m} = zeros (size(space.dvr{1}));
    end
    
    % Transform to adiabatic representation using Hermitian conjugate of eigenvector matrix 
    for m = 1:hamilt.coupling.n_eqs
        for n = 1:hamilt.coupling.n_eqs
            psi.dvr{m} = psi.dvr{m} + conj(hamilt.eig_vec{n,m}).*psi.dia{n};
        end
    end

    %% Adiabatic to diabatic ("back transformation")
elseif strcmpi(direction,'adi2dia')

    % Retrieve diabatic wavefunctions and potential matrices (cell arrays!)
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr{m}    = psi.dia{m};
        for n = m:hamilt.coupling.n_eqs
            if ~isempty (hamilt.pot{m,n}.dia)
                hamilt.pot{m,n}.dvr = hamilt.pot{m,n}.dia;
            else
                hamilt.pot{m,n}.dvr = [];
            end
        end
    end

end
