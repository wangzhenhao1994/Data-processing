%----------------------------------------------------------------------
%
% This function sets up the grids for the coordinates in position and
% "momentum" space (DVR and FBR, respectively) as well as the weights
% needed for integration in DVR.
%
% If input variable "state" is an object of class "wave", i. e.
% for fully quantum mechanical simulations, this function also
% constructs N-dim grids as direct products of the one-dim grids.
%
%----------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function init (state)
global space 

% Dimensionality of problem
space.n_dim = length(space.dof);

if isa(state,'ket') % this includes subclass 'rho'
    return
end

if isa(state,'wave') || ( isa(state,'traj') && space.n_dim<3 )
    
    % Total number of grid points
    space.n_tot = 1;
    for k = 1:space.n_dim
        space.n_tot = space.n_tot * space.dof{k}.n_pts;
    end
    
    % Some logging header
    prt.disp ( '***************************************************************' )
    prt.disp ( 'Setting up a direct product of the one-dimensional grids       ' )
    prt.disp ( '***************************************************************' )
    prt.disp ( ' ' )
    prt.disp ( [ 'Dimensionality of problem          : ' int2str(space.n_dim) ] )
    prt.disp ( [ 'Total number of grid points        : ' int2str(space.n_tot) ] )
    prt.disp ( ' ' )
    
    % Build one-dimensional grids in DVR and FBR space as cell arrays
    weight = cell(space.n_dim,1);
    x_grid = cell(space.n_dim,1);
    p_grid = cell(space.n_dim,1);
    
else
    
    % Some logging header
    prt.disp ( '***************************************************************' )
    prt.disp ( 'Setting up one-dimensional grids for each degree of freedom    ' )
    prt.disp ( 'Used for visualization only!                                   ' )
    prt.disp ( '***************************************************************' )
    prt.disp ( ' ' )
    
end

% Main loop over all degrees of freedom
for k = 1:space.n_dim
    
    % For Q/C simulations only: So far, only FFT grids are permitted
    if isa (state,'traj')
        if ~isa(space.dof{k},'dof.fft')
            prt.disp (['Found a non-FFT grid for dof : ' int2str(k)])
            prt.error ('This is (currently!) not allowed for (quantum-)classical dynamics!')
        end
    end
    
    % Build one-dimensional grids in DVR and FBR space  
    space.dof{k}.dof = k; % Tell each DVR which degree of freedom it represents.
    init_grid(space.dof{k});
    disp     (space.dof{k})
    
    % For Q/M simulations only: Save grids/weights
    if isa(state,'wave') || ( isa(state,'traj') && space.n_dim<3 ) 
        weight{k} = space.dof{k}.weight;
        x_grid{k} = space.dof{k}.x_grid;
        p_grid{k} = space.dof{k}.p_grid;
    end
end

% Remainder: For Q/M simulations only
if isa (state,'traj') && space.n_dim>2
    return
end

% For multidimensional calculations, construct multi-dimensional grids
if space.n_dim == 1 % Note that ndgrid(x) has the same effect as ndgrid(x,x)
    space.weight  = weight{1};
    space.dvr{1} = x_grid{1};
    space.fbr{1} = p_grid{1};
else
    [       tmpweights{1:space.n_dim}] = ndgrid(weight{1:space.n_dim});
    [space.dvr{1:space.n_dim}] = ndgrid(x_grid{1:space.n_dim});
    [space.fbr{1:space.n_dim}] = ndgrid(p_grid{1:space.n_dim});
end

% Finally, put all the weights together to get a matrix that has to be applied
% when summing over the position grid.
if space.n_dim > 1
    space.weight = tmpweights{1};
    for k = 2:space.n_dim
        space.weight = space.weight .* tmpweights{k};
    end
end

