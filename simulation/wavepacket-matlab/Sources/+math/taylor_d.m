%********************************************************************
%
% NEGATIVE gradient of Taylor series expansion (diagonal in N dimensions)
%
% df(R)   inf   N    c_jk            j-1            
% ----- = Sum  Sum  ------ ( R  - h )         
%  d R    j=1  k=1  (j-1)!    k    k            
%
%********************************************************************

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2016 Burkhard Schmidt
%
% see the README file for license details.

function df = taylor_d (R,hshift,coeffs)

% Preallocate gradients by setting them to zero
df = cell(length(hshift),1);
for k = 1:length(hshift)
    df{k} = zeros(size(R{1}));
end
    
% Return if no Taylor series coefficients available
if isempty(coeffs)
    return
end

% Size of coefficient matrix
[nrow,ncol]=size(coeffs);
if ncol~=length(hshift)
    prt.error ('Wrong length of Taylor coefficient vectors')
end

% Loop over expansion orders
for j=1:nrow
    fj = factorial(j);
    
    % Summing up contributions from each component of position vector
    for k = 1:ncol
        df{k} = df{k} - j * (R{k}-hshift(k)).^(j-1) * coeffs(j,k) / fj;
    end
end