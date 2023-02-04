%--------------------------------------------------------------------------
%
% Calculate matrix elements ("sandwiches") <bra|operator|ket> 
% by DVR (quadrature) methods
%
% where "bra" and "ket" are objects of class "wave" and
% where "operator" is an object of a suitable class, e.g. pot, dip, ...
%
%--------------------------------------------------------------------------


% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017 Burkhard Schmidt
%               2011 Ulf Lorenz, Boris Schaefer-Bung, Burkhard Schmidt
%               2012 Jeremy Rodriguez, Burkhard Schmidt, Ulf Lorenz
%
% see the README file for license details.

function retval = sandwich (bra, operator, ket)
global space hamilt

retval = 0;

% If operator is a cell array (for coupled channels)
if iscell (operator)
    
    for m = 1:hamilt.coupling.n_eqs
        % diagonal couplings count only once
        if ~isempty(operator{m,m}.dvr)
            retval = retval + sum(conj(bra{m}(:)) .* operator{m,m}.dvr(:) .* ket{m}(:) .* space.weight(:));
        end
        
        % the offdiagonal coupling maps both off-diagonal elements on one term only
        for n = m+1:hamilt.coupling.n_eqs
            if ~isempty(operator{m,n}.dvr)
                retval = retval ...
                    + sum(conj(bra{m}(:)) .*      operator{m,n}.dvr(:)  .* ket{n}(:) .* space.weight(:)) ...
                    + sum(conj(bra{n}(:)) .* conj(operator{m,n}.dvr(:)) .* ket{m}(:) .* space.weight(:));
            end
        end
    end
    
% If operator acts in the same way for each of the (uncoupled) channels)
else
    
    for m = 1:hamilt.coupling.n_eqs
        retval = retval + sum(conj(bra{m}(:)) .* operator.dvr(:) .* ket{m}(:) .* space.weight(:));
    end

end

end

