%--------------------------------------------------------------------------
%
% Solve eigenproblem for the instantaneous (effective) potential matrix
%
% For the special case of two coupled Schrödinger equations, 
% the eigenvalue problem can be solved analytically 
% 
% Diabatic: V = ( alpha  beta )
%               ( beta  gamma )
%
% delta = (alpha-gamma)/2
%   eta = (alpha+gamma)/2
%   rho = sqrt(delta^2+beta^2)>0
% theta = atan ( beta/delta )
%
% Adiabatic: E = (eta-rho    0    )
%                (   0    eta+rho )
%
% dia2adi: S = ( -sin(theta/2)  +cos(theta/2) )
%              ( +cos(theta/2)  +sin(theta/2) )
%
% The first and second column of S contain the eigenvectors
% corresponding to the lower (label 1) and upper (label 2)
% eigenvalue contained in the diagonal entries of matrix E.
% Note that in this case we have S^-1 = S^+ = S, i.e. the
% 'dia2adi' and 'adi2dia' transformation are described by 
% the same matrix.
%
% Obviously, this transformation results in double-valued
% adiabatic states upon increasing theta from 0 to 2*pi, e.g.
% upon encircling a conical intersection of two adiabatic
% potential energy surfaces. This can be remedied by intro-
% ducing complex-valued adiabatic states, i.e. by multiplying
% S with exp(i*phi) or multiplying S^+ with exp(-i*phi). The
% geometric phase has to be chosen as phi = n*theta/2 where
% n is an odd integer. Note that this is equivalent to the 
% introduction of an appropriately chosen vector potential 
% governing the dynamics of wavepackets along the adiabatic 
% potential energy surfaces. (Alternatively, one could compensate
% the double-valuedness of the adiabatic states by enforcing 
% doubled-valued boundary conditions for the wavefunctions, too.)
%
% Input argument "e" is the electric field (vector)
%
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008,2010 Ulf Lorenz
%
% see the README file for license details.

function pot_eig_2 ( e, save_trafo )
global hamilt space

% Bare potential: state 1
if ~isempty ( hamilt.pot{1,1}.dvr ) 
    alpha = hamilt.pot{1,1}.dvr;
else % free particle
    alpha = zeros ( size ( space.dvr{1} ) );
end

% Bare potential: state 2
if ~isempty ( hamilt.pot{2,2}.dvr ) 
    gamma = hamilt.pot{2,2}.dvr;
else % free particle
    gamma = zeros ( size ( space.dvr{1} ) );
end

% Permanent dipole moments along p = x|y
if isfield(hamilt, 'dip')
    for p = 1:length(hamilt.dip)
        if abs(e(p))>0
            if ~isempty (hamilt.dip{p})
                if ~isempty ( hamilt.dip{p}{1,1}.dvr )
                    alpha = alpha - e(p) * hamilt.dip{p}{1,1}.dvr;
                end
                if ~isempty ( hamilt.dip{p}{2,2}.dvr )
                    gamma = gamma - e(p) * hamilt.dip{p}{2,2}.dvr;
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
                            alpha = alpha - e(p)*e(q)/2 * hamilt.pol{p,q}{1,1}.dvr;
                        end
                        if ~isempty ( hamilt.pol{p,q}{2,2}.dvr  )
                            gamma = gamma - e(p)*e(q)/2 * hamilt.pol{p,q}{2,2}.dvr ;
                        end
                    end
                end
            end
        end
    end
end

% Three possible coupling mechanisms
beta  =  zeros (size(space.dvr{1}));

% Diabatic potential coupling
if ~isempty ( hamilt.pot{1,2}.dvr ) 
    beta  =  beta + hamilt.pot{1,2}.dvr;
end

% Transition dipole moments along p = x|y
if isfield(hamilt, 'dip')
    for p = 1:length(hamilt.dip)
        if abs(e(p))>0
            if ~isempty (hamilt.dip{p})
                if  ~isempty(hamilt.dip{p}{1,2}.dvr)
                    beta = beta - e(p) * hamilt.dip{p}{1,2}.dvr;
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
                    if ~isempty (hamilt.pol{p,q}.dvr)
                        if  ~isempty(hamilt.pol{p,q}{1,2}.dvr)
                            beta = beta - e(p)*e(q)/2 * hamilt.pol{p,q}{1,2}.dvr;
                        end
                    end
                end
            end
        end
    end
end

% Sum, difference of diabatic potentials
eta = (alpha + gamma)/2; % half trace
dlt = (alpha - gamma)/2; % half gap

% Polar coordinates: "Radius"
rho = sqrt ( dlt.^2 + beta.^2 );

% Adiabatic potential energy surfaces: 1=lower, 2=upper
hamilt.eig_val    = cell (2,1); % cell vector
hamilt.eig_val{1} = eta - rho;
hamilt.eig_val{2} = eta + rho;

% Adiabatic<=>diabatic transformation matrices
% Columns of eig_vec are the right eigenvectors corresponding to eig_val
if save_trafo 
    theta = atan2 ( beta , dlt ); % mixing angle
    hamilt.eig_vec    = cell (2); % cell matrix
    hamilt.eig_vec{1,1} = - sin(theta/2) .* exp(1i*theta/2);
    hamilt.eig_vec{1,2} = + cos(theta/2) .* exp(1i*theta/2);
    hamilt.eig_vec{2,1} = + cos(theta/2) .* exp(1i*theta/2);
    hamilt.eig_vec{2,2} = + sin(theta/2) .* exp(1i*theta/2);
end
end

