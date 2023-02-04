%--------------------------------------------------------------------------
%
% Propagate quantum state vectors attached to trajectories
% The time evolution is given by the Schroedinger equation
% 
%    d           i   ^         
%   -- c(t) = - ---- H(t) c(t) 
%   dt -        hbar =    -    
% 
% for a time-step TAU using second (or higher) order differencing
% 
%                            i       ^            i    ^3
%   c(t+tau) = c(t-tau) - 2 ---- tau H c(t) + -------- H  c(t) - ...
%                           hbar     = -      3 hbar^3 =  -
% 
% This version for diabatic representation 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20.. Burkhard Schmidt's group
%
% see the README file for license details.

function tdse_dia (obj)
global hamilt time

% Initialize: Propagate state vectors psi backward using first order method
if isempty(obj.psi_old{1})
    for m=1:hamilt.coupling.n_eqs
        obj.psi_old{m} = obj.psi{m};
        for n=1:hamilt.coupling.n_eqs
            obj.psi_old{m} = obj.psi_old{m} + 1i * obj.ham{m,n} .* obj.psi{n} * time.steps.s_delta;
        end
    end
    
else
    
    H       = zeros ( obj.n_p , hamilt.coupling.n_eqs , hamilt.coupling.n_eqs );
    H_trans = zeros ( size(H) ); % H^T
    for m=1:hamilt.coupling.n_eqs
        for n=1:hamilt.coupling.n_eqs
            H(:,m,n)        = -1i * time.steps.s_delta * obj.ham{m,n};
            H_trans(:,n,m)  = H(:,m,n);
        end
    end
    
    % Has to be odd
    maxdegree = 11;
    
    if(mod(maxdegree,2)==0)
        prt.error ( 'maxdegree has to be odd' )
    end
    
    H_power             = cell(maxdegree,1);    % the j-th cell entry will be H^j
    H_power_trans       = cell(size(H_power));  % the j-th cell entry will be (H^j)^T
    H_power{1}          = H;
    H_power_trans{1}    = H_trans;
    
    % Computation of H^2
    if(maxdegree > 2)
        H_power{2}          = zeros ( size(H) );
        H_power_trans{2}    = zeros ( size(H) );
        for m=1:hamilt.coupling.n_eqs
            for n=1:hamilt.coupling.n_eqs
                % Computation of H^2 and (H^2)^T
                % A * B = element-wise product of A^T and B with 
                % the summmation dimension over the right dimension
                H_power{2}(:,m,n)       = sum( H_trans(:,:,m) .* H(:,:,n) , 2 );
                H_power_trans{2}(:,n,m) = H_power{2}(:,m,n);
            end
        end
    end
    
    % Computation of H^j for odd j
    for j = 3:2:maxdegree
        H_power{j}          = zeros ( size(H) );
        H_power_trans{j}    = zeros ( size(H) );
        for m=1:hamilt.coupling.n_eqs
            for n=1:hamilt.coupling.n_eqs
                H_power{j}(:,m,n)       = sum( H_power_trans{j-2}(:,:,m) .* H_power{2}(:,:,n) , 2 );
                H_power_trans{j}(:,n,m) = H_power{j}(:,m,n);
            end
        end
    end
    
    % Computation of the exponential series
    for m=1:hamilt.coupling.n_eqs
        obj.psi_new{m} = obj.psi_old{m};
        for n=1:hamilt.coupling.n_eqs
            for j = 1:2:maxdegree
                obj.psi_new{m} = obj.psi_new{m} + 2 / factorial(j) * H_power{j}(:,m,n) .* obj.psi{n};
            end
        end
    end

    % Save state vectors psi; get ready for next step
    for m=1:hamilt.coupling.n_eqs
        obj.psi_old{m} = obj.psi    {m};
        obj.psi    {m} = obj.psi_new{m};
    end
    
end
end
