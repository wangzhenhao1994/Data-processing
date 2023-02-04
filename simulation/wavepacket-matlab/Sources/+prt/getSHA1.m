%--------------------------------------------------------------------------
% Get SHA-1 hash value (used within GIT version control) for a given path 
%
% Using GIT command line. If not installed, a "[]" will be returned.
%
% By passing the option "--abbrev-commit" to the git log command, 
% the output will use shorter values but keep them unique; 
% it defaults to using seven characters but makes them longer 
% if necessary to keep the SHA-1 unambiguous (from the GIT bible)
%
% By passing the option "--pretty=oneline" to the git log command, 
% the output (SHA-1 hash value and commit message) will be in one line. 
%
% By passing the argument "-1" to the git log command, 
% only the latest commit will be displayed.
%
%--------------------------------------------------------------------------
function retval = getSHA1 (pathstr)

command = ['git log --abbrev-commit --pretty=oneline -1 ' pathstr];
[status,cmdout] = system (command);
if status ~= 0
    retval = [];
    return
else
    retval = [cmdout(1:7) '...'];
end