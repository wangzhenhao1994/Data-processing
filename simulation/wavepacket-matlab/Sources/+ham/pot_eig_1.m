%--------------------------------------------------------------------------
%
% Get instantaneous (effective) potential energy curve
% for a single channel Schrödinger equation
%
% Input argument "e" is the electric field (vector)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008,2010 Ulf Lorenz
%
% see the README file for license details.

function pot_eig_1 ( e, save_trafo )
global hamilt space

% Bare potential
if ~isempty ( hamilt.pot{1,1}.dvr ) 
    dressed = hamilt.pot{1,1}.dvr;
else % free particle
    dressed = zeros (size ( space.dvr{1} ) );
end
    
% Dipole moments along p = x|y
if isfield(hamilt, 'dip')
    for p = 1:length(hamilt.dip)
        if abs(e(p))>0
            if ~isempty (hamilt.dip{p})
                if ~isempty ( hamilt.dip{p}{1,1}.dvr )
                    dressed = dressed - e(p) * hamilt.dip{p}{1,1}.dvr;
                end
            end
        end
    end
end

% Polarizabilities along p = x|y
if isfield(hamilt, 'pol')
    for p = 1:size(hamilt.pol,1)
        if abs(e(p))>0
            for q = p:size(hamilt.pol,2)
                if abs(e(q))>0
                    if ~isempty (hamilt.pol{p,q})
                        if ~isempty ( hamilt.pol{p,q}{1,1}.dvr )
                            dressed = dressed - e(p)*e(q)/2 * hamilt.pol{p,q}{1,1}.dvr;
                        end
                    end
                end
            end
        end
    end
end

% Eigenvalues (effective potential)
hamilt.eig_val    = cell (1);
hamilt.eig_val{1} = dressed;

% Transformation matrix: trivial
if save_trafo 
    hamilt.eig_vec    = cell (1);
    hamilt.eig_vec{1,1} = ones(size(dressed));
end

end