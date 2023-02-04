%--------------------------------------------------------------------------
%
% General information: Names/versions of program, user, host, etc
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-20xy Burkhard Schmidt
%               2007-2014 Ulf Lorenz
%
% see the README file for license details.

function init (filename)
global info

% Path and file name
[pathstr, filestr, ~] = fileparts (filename);

% Program name and version number
info.package = 'WavePacket';
info.program = filestr;

% Get SHA-1 hash value: WavePacket installation
workdir = pwd;
cd (pathstr)
sha = prt.getSHA1 ('./');
if ~isempty(sha)
    info.version1 = sha;
else
    info.version1 = 'V7.0.0';
end
cd (workdir)

% Get SHA-1 hash value: Current working directory
sha = prt.getSHA1 ('./');
if ~isempty(sha)
    info.version2 = sha;
else
    info.version2 = '(unversioned)';
end

% Set/unset these two lines before/after each release of a new version
info.version1 = 'V7.0.0';
info.version2 = 'V7.0.0';

% Get MATLAB or Octave release (date)
if strcmpi(info.system,'Matlab')
    info.release = version ('-release');
elseif strcmpi(info.system,'Octave')
    info.release = version ('-date');
end

% Get MATLAB or Octave user name 
if strcmpi(info.system,'Matlab')
    result = license('inuse', 'matlab');
elseif strcmpi(info.system,'Octave')
    result = license('inuse', 'octave');
end
info.user_name = result.user;


% Get host name
if isunix 
    [status,result] = unix('hostname');
    if ~status
        info.host_name = strcat(result,' (Unix)');
    else
        info.host_name = 'Running on Unix';
    end
elseif ispc
    [status,result] = dos('hostname');
    if ~status
        info.host_name = strcat(result,' (Windows)');
    else
        info.host_name = 'Running on Windows';
    end
end
    
% Get path and file name 
info.path_name = pwd;
addpath(info.path_name);

% Open log file (or create a new one) for writing in text mode
% Discard existing contents, if any.
filename = fullfile(info.path_name, strcat(info.program, '.log'));
fclose('all'); % needed when qm_xyz is called inside a loop
info.stdout = fopen(filename, 'wt');

if info.stdout == -1
    error(strcat('Could not open log file ',filename));
end

% Output program name/version etc
prt.disp ( ' ' )
prt.disp ( '***************************************************************')
prt.disp ( 'About this WavePacket simulation' )
prt.disp ( '***************************************************************')
prt.disp ( ' ' )
prt.disp (['Program package : ', info.package] )
prt.disp (['Program name    : ', info.program] )
prt.disp (['SHA1 hash value : ', num2str(info.version1),' (source codes)'])
prt.disp (['SHA1 hash value : ', num2str(info.version2),' (working dir.)'])
prt.disp ( ' ' )
prt.disp (['Path name       : ', info.path_name] )
prt.disp (['User name       : ', info.user_name] )
prt.disp (['Host name       : ', info.host_name] )
if strcmpi(info.system,'Matlab')
  prt.disp (['MATLAB release  : ', info.release  ] )
elseif strcmpi(info.system,'Octave')
  prt.disp (['Octave release  : ', info.release  ] )
end
prt.disp ( ' ' )

% Initialize stopwatch timer; Output clock/date/time
info.start_time = cputime;
prt.clock;

% Output GPL and copyright info
prt.disp ( '***************************************************************')
prt.disp ( 'Licensing information ' )
prt.disp ( '***************************************************************')
prt.disp ( ' ' )
prt.disp ( 'This program is subject to the GNU General ')
prt.disp ( 'Public License v2, see http://www.gnu.org or the ')
prt.disp ( 'root path for the license text')
prt.disp ( ' ' )
prt.disp ( ' (C) 2004-2019 Burkhard Schmidt, Ulf Lorenz, FU Berlin' )
prt.disp ( '     http://sourceforge.net/projects/matlab.wavepacket.p' )
prt.disp ( ' ' )
