%-------------------------------------------------------------------------
%
% Solve the time-independent Schroedinger equation to get eigenstates
% and energies in pseudospectral representation using DVR/FBR techniques 
%
% Part 2/3: Optionally enforce symmetry restrictions. 
%           So far, this is restricted to a 1-D grid.

%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function symm (obj)

global hamilt space

% Symmetry adaption is optional
if isempty(obj.symmetry) || strcmpi(obj.symmetry,'n')
    return
else
    if space.n_dim > 1
        prt.error ('So far, symmetry adaption has been implemented only for 1-dim systems')
    end
end

% Total number (N) of grid points should be even!
N = space.n_tot;
if mod(N,2)
    prt.error ('Symmetry adaption only for even number of grid points')
end

% Mapping the full problem to a symmetry-adapted basis
switch class(space.dof{1})
    
    case {'dof.legendre','dof.hermite'}
        main = eye ( N/2);
        anti = fliplr(main);
        switch lower(obj.symmetry)
            case 'g' % N/2 rows in transformation matrix
                prt.disp ('Even (g) parity (A_g irrep of C_i group) - GAU-LEG-HER')
                obj.transform = [main anti];
            case 'u' % N/2 rows in transformation matrix
                prt.disp ('Even (g) parity (A_g irrep of C_i group) - GAU-LEG-HER')
                obj.transform = [main -anti];
            otherwise
                prt.error ('Symmetry should be "g" (=even) or "u" (=odd) for GAU-LEG-HER')
        end
        
    case 'dof.fft' % Note absence of grid point at space.dof{1}.r_max
        switch obj.symmetry
            case 'g' % N/2+1 rows in transformation matrix
                prt.disp ('Even (g) parity (A_g irrep of C_i group) - FFT')
                main = eye ( N/2-1 );
                anti = fliplr(main);
                init = zeros(1,N); init(1)=1;
                last = zeros(1,N); last(N/2+1)=1;
                fill = zeros(N/2-1,1);
                obj.transform = [init; fill main fill anti; last];
            case 'u' % N/2-1 rows in transformation matrix
                prt.disp ('Odd (u) parity (A_u irrep of C_i group) - FFT')
                main = eye ( N/2-1 );
                anti = fliplr(main);
                fill = zeros(N/2-1,1);
                obj.transform = [fill main fill -anti];
            case 'A1'
                prt.disp ('1st (A_1) irrep of C_2v group - FFT')
                main = eye ( N/4-1 );
                anti = fliplr(main);
                init = zeros(1,N); init(1)=1; init(N/2+1)=1;
                last = zeros(1,N); last(N/4+1)=1;last(3*N/4+1)=1;
                fill = zeros(N/4-1,1);
                obj.transform = [init; fill main fill anti fill main fill anti; last];
            case 'B1'
                prt.disp ('2nd irrep of C_2v group - FFT')
                main = eye ( N/4-1 );
                anti = fliplr(main);
                init = zeros(1,N); init(1)=1; init(N/2+1)=-1;
                fill = zeros(N/4-1,1);
                obj.transform = [init; fill main fill -anti fill -main fill anti];
            case 'B2'
                prt.disp ('3rd irrep of C_2v group - FFT')
                main = eye ( N/4-1 );
                anti = fliplr(main);
                last = zeros(1,N); last(N/4+1)=1;last(3*N/4+1)=-1;
                fill = zeros(N/4-1,1);
                obj.transform = [fill main fill anti fill -main fill -anti; last];
            case 'A2'
                prt.disp ('4th irrep of C_2v group - FFT')
                main = eye ( N/4-1 );
                anti = fliplr(main);
                fill = zeros(N/4-1,1);
                obj.transform = [fill main fill -anti fill main fill -anti];
            otherwise
                prt.error ('Symmetry should be "g" (=even) or "u" (=odd) or A1,A2,B1,B2 for FFT')
        end
        
end
    
% normalize rows of transformation matrices
for j=1:size(obj.transform,1)
    obj.transform(j,:) = obj.transform(j,:) ...
        / norm(obj.transform(j,:));
end

% for more than one surface, we just create a block-diagonal matrix, where each
% block looks the same, and transforms the wave function on the corresponding
% surface to the reduces representation.
if hamilt.coupling.n_eqs > 1
    newtrafo = zeros(space.n_tot/2 + hamilt.coupling.n_eqs, space.n_tot);
    for m = 1:hamilt.coupling.n_eqs
        nrows = size(obj.transform,1);
        newtrafo((m-1)*nrows+1:m*nrows, (m-1)*N+1:m*N) = obj.transform;
    end
    obj.transform = newtrafo;
end
    
    % transform Hamiltonian matrix
    obj.matrix = obj.transform*obj.matrix*obj.transform';
end

