%--------------------------------------------------------------------------
%
% Shorten a string, by taking only the last and first 21 characters
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2014 - 2017 Burkhard Schmidt
%               2011 Ulf Lorenz, Boris Schaefer-Bung, Burkhard Schmidt
%               2012 Jeremy Rodriguez, Burkhard Schmidt, Ulf Lorenz
%               
%
% see the README file for license details.

function str_output = shorten( str_input )

L = length(str_input);

if L>45
    str1 = str_input(1:21);
    str2 = str_input(L-20:L);
    
    str_output = [str1,' ... ',str2];
    
else
    
    str_output = str_input;

end

