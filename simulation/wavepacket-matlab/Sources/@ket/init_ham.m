%--------------------------------------------------------------------------
%
% Read A, B, N, C, D matrices from data file: usually ket_0.mat
% for use with bilinear control systems:
% A   : system
% B, N: control
% C, D: observe
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-20.. Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2011 Boris Schaefer-Bung, Ulf Lorenz
%
% see the README file for license details.

function init_ham (state)

load_0 ( state, 1 )

