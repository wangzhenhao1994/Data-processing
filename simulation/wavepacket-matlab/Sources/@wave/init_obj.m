%------------------------------------------------------------------------------
%
% Set up the initial state of a wave function object 
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function init_obj (obj)
global hamilt space time

% Create cell vectors (column vectors)
obj.dvr = cell (hamilt.coupling.n_eqs,1);
obj.dia = cell (hamilt.coupling.n_eqs,1);
obj.adi = cell (hamilt.coupling.n_eqs,1);
obj.fbr = cell (hamilt.coupling.n_eqs,1);
obj.ini = cell (hamilt.coupling.n_eqs,1);
obj.new = cell (hamilt.coupling.n_eqs,1);
obj.old = cell (hamilt.coupling.n_eqs,1);
obj.sum = cell (hamilt.coupling.n_eqs,1);
obj.wig = cell (hamilt.coupling.n_eqs,1);

% Initialize the wave functions with ZEROs
for m = 1:hamilt.coupling.n_eqs
	obj.dvr{m} = zeros(size(space.weight));
end

% Specify the initial state EITHER as correlated OR as product state
if isfield(time, 'corr') && isfield(time, 'dof')
    prt.error ('Specify the initial state EITHER as correlated OR as product state')
end

% Fully correlated initial wave functions
if isfield(time, 'corr')
    if ~isempty(time.corr)
        prt.disp ('***************************************************************')
        prt.disp ('Fully correlated initial wavefunction  ')
        init ( time.corr );
        disp ( time.corr ); prt.disp ( ' ' )
        wave ( time.corr , obj);
    else
        prt.error ('No specification for correlated initial wavefunction found')
    end
end

% Direct poduct of 1D functions
if isfield(time, 'dof')
    if space.n_dim>1
        prt.disp ('***************************************************************')
        prt.disp ('Initial wavefunction as direct product of 1-dim states         ')
        prt.disp ('***************************************************************')
        prt.disp (' ')
    end
    
    if length(time.dof) ~= space.n_dim
        prt.error ('Wrong dimensionality of initial product state')
    end
    
    obj.dvr{1} = ones(size(space.dvr{1}));
    for k = 1:space.n_dim
        prt.disp ('***************************************************************')
        if space.n_dim>1
            prt.disp (['Initial state along degree of freedom : ' space.dof{k}.label])
        end
        if ~isempty(time.dof{k})
            time.dof{k}.dof = k; % Tell each d-o-f about its index
            init ( time.dof{k} );
            disp ( time.dof{k} ); prt.disp ( ' ' )
            wave ( time.dof{k} );
            obj.dvr{1} = obj.dvr{1} .* time.dof{k}.dvr;
        else
            prt.error ('No specification for this d-o-f found')
        end
    end
    
end

%% Normalize the wavefunction
if hamilt.coupling.ini_norm
    obj.dvr = wave.normalize(obj.dvr);
end

%% Population of coupled states
if hamilt.coupling.n_eqs>1
    
    % Channels initially populated according to individual coefficients
    if ~isempty(hamilt.coupling.ini_coeffs)
        
        % Check number of coefficients
        if length(hamilt.coupling.ini_coeffs) ~= hamilt.coupling.n_eqs
            prt.error('Wrong number of initial coefficients')
        end

        % Diabatic initial states
        for m=hamilt.coupling.n_eqs:-1:1
            obj.dvr {m} = obj.dvr{1} * hamilt.coupling.ini_coeffs(m);
        end

        % If required, transform initial state from adiabatic to diabatic representation
        if strcmpi(hamilt.coupling.ini_rep,'adi')
            diabatic (obj);
        end

    else
        prt.disp ('***************************************************************')
        prt.disp ('Initial populations taken from initial function.  ')
        prt.disp ('***************************************************************')
        prt.disp (' ')

        % Has to be separated, otherwise isempty(...) gives an error
        for m = 1:min(hamilt.coupling.n_eqs, numel(obj.dvr))
            if isempty(obj.dvr{m})
                obj.dvr{m} = zeros(size(space.dvr{1}));
            end
        end
        for m = numel(obj.dvr)+1:hamilt.coupling.n_eqs
            obj.dvr{m} = zeros(size(space.dvr{1}));
        end
        
    end
    
end


%% Save initial wavefunction for later use (->autocorrelation)
for m = 1:hamilt.coupling.n_eqs
    obj.ini{m} = obj.dvr{m};
end
