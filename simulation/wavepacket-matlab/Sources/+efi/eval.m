%--------------------------------------------------------------------------
%
% External electric field as a sequence of shaped pulses with envelope f
% and constant or (linearly or quadratically) chirped carrier frequencies 
% with user-specified polarizations in the x/y plane and
% possibly with phase shifts.
%
%            N
% E(time) = sum f (time) * cos ( omega  * time + phase )
%           i=1  i                    i               i
%
% Returning a cell vector of length two. Each entry contains a vector of 
% field values for any given vector (timesteps) of time values.
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function efield = eval (timesteps)
global time

% Preallocate
efield{1} = zeros(size(timesteps));
efield{2} = zeros(size(timesteps));

% Loop over pulses
for p = 1:length(time.pulse)
    if ~isempty(time.pulse{p})
        
        % Get envelopes and oscillations via methods of the pulse objects
        shape = envelope  (time.pulse{p},timesteps);
        vibra = oscillate (time.pulse{p},timesteps);
        
        % Find those timesteps where pulse is "on"
        where = find( abs(timesteps-time.pulse{p}.delay) <= time.pulse{p}.length/2 );
        
        % Dressed state picture (Floquet)
        if time.efield.dressed % Return half envelope
            for xy=1:2
                efield{xy}(where) = efield{xy}(where) + ...
                    0.5 * shape(where) * time.pulse{p}.strength(xy);
            end
            
            % Bare state picture
        else % Return oscillating field
            for xy=1:2
                efield{xy}(where) = efield{xy}(where) + ...
                    shape(where) .* vibra(where) * time.pulse{p}.strength(xy);
            end
        end
        
    end
    
end
end
