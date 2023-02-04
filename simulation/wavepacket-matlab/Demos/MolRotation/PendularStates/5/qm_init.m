% Copyright (C) 2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (~)
global hamilt plots space

prt.disp ( '***************************************' )
prt.disp ( 'Pendular States for \Delta \omega = 80 ' )
prt.disp ( '***************************************' )

% Spatial discretization
space.dof{1}       = dof.legendre;       % Legendre polynomials in cos theta
space.dof{1}.label = 'cos \Theta';
space.dof{1}.R_0 = 1;                    % constant value for R
space.dof{1}.m_0 = 0;                    % minor quantum number, fixed to 0
space.dof{1}.l_max = 201;                % maximum angular momentum/ number of points
space.dof{1}.mass = 0.5;                 % adjusted mass

% Orientation: cosine functiom
hamilt.amo{1} = amo.cosine;
hamilt.amo{1}.exp = 1;

% Alignment: cosine^2 function
hamilt.amo{2} = amo.cosine;
hamilt.amo{2}.exp = 2;

% Hamiltonian operator 
hamilt.pot{1,1}        = pot.taylor;     % Taylor series in cos theta
hamilt.pot{1,1}.coeffs = [0;-2*80];      % Force constant

% Select eigen/values/functions
hamilt.eigen.start     =  0;             % Lower index
hamilt.eigen.stop      = 15;             % Upper index
hamilt.eigen.symmetry  = 'n';            % may be used if potential symmetric

% Plot densities
plots.density         = vis.polar;       % polar plot

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 300;                % Set maximum for energy plot
