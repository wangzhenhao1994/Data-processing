%--------------------------------------------------------------------------
%
% Check stability of (sparse) system matrix A: 
% real part of all eigenvalues should be negative
%
%--------------------------------------------------------------------------

function check_stable ( obj, label )

% Construct matrix and calculate eigenvalues
matrix = obj.A; % + obj.A';
if nargin>2
    for d=1:length(obj.N)
        matrix = matrix + obj.N{d} * obj.N{d}';
    end
end
eigval=sort(real(eig(full(matrix))),'descend');

% Check for non-negative real part
epsilon = eps('single');
nnrp = nnz(eigval>=-epsilon);
prt.disp([num2str(nnrp) ' non-negative eigenvalues found for ' label ])

% Display eigenvalue with largest real part
if nnrp
    if nnrp < 10
        prt.disp('Eigenvalue(s) with largest real part(s)')
        for k=1:nnrp
            prt.disp ([int2str(k) ' : ' num2str(eigval(k))])
        end
    else
        for k=1:5
            prt.disp ([int2str(k) ' : ' num2str(eigval(k))])
        end
        prt.disp('... ...')
        for k=nnrp-4:nnrp
            prt.disp ([int2str(k) ' : ' num2str(eigval(k))])
        end
        
    end
end

end

