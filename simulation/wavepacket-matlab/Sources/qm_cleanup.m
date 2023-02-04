%------------------------------------------------------------------------------
%
% Cleans up the calculation things (mainly closes log file for now).
%
% https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_cleanup
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt
%               2011 Ulf Lorenz
%
% see the README file for license details.

function qm_cleanup()

% For relaxation with 'cheby_imag': save bound state
% To be used in next relaxation (projecting out)
global hamilt info time state
if isfield (time,'propa')
    if isa (time.propa,'tmp.wave.cheby_imag')
        
        % Number of columns of cell array state.bound
        n_col = size (state.bound,2);
        
        % Adding a column to cell array state.bound
        for m=1:hamilt.coupling.n_eqs
            state.bound{m,n_col+1} = state.dvr{m};
        end
    end
end
% Close log file
fclose (info.stdout);
info.stdout = [];

% Remove path
rmpath(info.path_name);
