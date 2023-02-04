%********************************************************************
%
% Taylor series expansion (diagonal in N dimensions)
%
%         inf   N   c_jk           j             
% f (R) = Sum  Sum  ---- ( R  - h )  + v        
%         j=1  k=1   j!     k    k            
%
%********************************************************************

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2016 Burkhard Schmidt
%
% see the README file for license details.

function fr = taylor (R,hshift,vshift,coeffs,verbose)

% Use constant offset to initialize function
fr = vshift * ones(size(R{1}));
    
% Confirm input variables
if verbose
    prt.disp ( [ 'Vertical shift (v)        : ' num2str(vshift) ] )
    prt.disp ( [ 'Horizontal shift (h)      : ' num2str(hshift) ] )
    prt.disp (' ')
end

% Return if no Taylor series coefficients available
if isempty(coeffs)
    return
end

% Size of coefficient matrix
[nrow,ncol]=size(coeffs);
if ncol~=length(R)
    prt.error ('Wrong length of Taylor coefficient vectors')
end

if length(hshift)~=length(R)
    prt.error ('Wrong length of horizontal shift vector')
end

if verbose
    prt.disp (['Taylor expansion order    : ' int2str(nrow)])
end

% Loop over expansion orders
for j=1:nrow
    fj = factorial(j);
    
    % Tabulate expansion coefficients
    if verbose
        prt.disp ( [ '==> ' int2str(j) '-th order coefficient(s) : ' num2str(coeffs(j,:)) ] )
    end
    
    % Summing up contributions from each component of position vector
    for k = 1:ncol
        fr = fr + (R{k}-hshift(k)).^j * coeffs(j,k) / fj;
    end
end

if verbose
    prt.disp (' ')
end