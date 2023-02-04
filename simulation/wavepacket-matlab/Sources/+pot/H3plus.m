% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%
% see the README file for license details.

classdef H3plus < pot.generic & handle
    
    properties (Access = public)
        
        c_dof       % angular coordinate
        
    end
    
    properties (Access = private)
        
        beta
        R_e
        Vmat
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = H3plus
            
            obj.beta = 1.3;
            obj.R_e  = 1.6501;
            
            obj.Vmat = zeros(8, 4, 3);
            
            % V_nm0
            obj.Vmat(:,:, 1) = [ 0      266219  44851  15068;
                130   -241851 -11820  37493;
                204603 131115  120688 0;
                -49832 -50919   23840  0;
                25002  50424   0      0;
                -2115   887     0      0;
                4346   0       0      0;
                -277    0       0      0];
            
            % V_nm1
            obj.Vmat(:,:,2) = [-6490 -3185   7605 0;
                88648 73273  0    0;
                -28688 104361 0    0;
                57028 0      0    0;
                9333  0      0    0;
                0     0      0    0;
                0     0      0    0;
                0     0      0    0];
            
            % V_nm2
            obj.Vmat(:, 1, 3) = [-339; -3238; 0; 0; 0; 0; 0; 0];
            
            obj.Vmat = obj.Vmat * 1e-6;
            
        end
        
        % Initialize potential: Set/check parameters
        function init (obj)
            global hamilt space
            if hamilt.coupling.n_eqs ~= 1
                prt.error ('This potential only for single Schrödinger equation')
            end
            
            if space.n_dim ~= 3
                prt.error ('This potential is for 3 degrees of freedom only')
            end
            
            if ~isa (space.dof{obj.c_dof}, 'dof.legendre')
                prt.error ('The potential requires Jacobi coordinates.')
            end
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            prt.disp ('***************************************************************')
            prt.disp ('High quality MRCI calculation for potential energy')
            prt.disp ('surface of  H3+                                   ')
            prt.disp ('See J. Chem. Phys. 84:891                         ')
            prt.disp ('                                                  ')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp ( [ 'angular coordinate : ' num2str(obj.c_dof) ] )
            prt.disp (' ')
        end

        % Evaluate  potential energy functions
        function V = V(obj,r)
            global space 
            % The transformation of the variables goes in three steps.
            % 1. We transform from Jacobi to Bond length coordinates
            % 2. The bond length coordinates are transformed in an exponential form
            % 3. Using these exponentiated coordinates, one can define the S3
            %    symmetry-adapted deformation coordinates
            %
            % The potential is then given as a power series in the deformation
            % coordinates.

            % First, define the angle and the Jacobi lengths
            if obj.c_dof == 1
                ang = r{1};
                j1  = r{2};
                j2  = r{3};
            elseif obj.c_dof == 2
                ang = r{2};
                j1  = r{1};
                j2  = r{3};
            else
                ang = r{3};
                j1  = r{1};
                j2  = r{2};
            end
            
            % Now go over to bond length coordinates
            r12 = j1;
            r13 = sqrt( j1.^2/4 + j2.^2 + 2 * j1/2 .* j2 .* ang );
            r23 = sqrt( j1.^2/4 + j2.^2 - 2 * j1/2 .* j2 .* ang );
            
            % Transform to the exponential form
            q12 = 1/obj.beta * ( 1 - exp(-obj.beta * (r12/obj.R_e-1)) );
            q13 = 1/obj.beta * ( 1 - exp(-obj.beta * (r13/obj.R_e-1)) );
            q23 = 1/obj.beta * ( 1 - exp(-obj.beta * (r23/obj.R_e-1)) );
            
            % Finally, transform to deformation coordinates
            sa = (q12 + q13 + q23) / sqrt(3);
            sx = (2*q12 - q23 - q13) / sqrt(6);
            sy = (q23 - q13) / sqrt(2);
            
            % The coordinates x/y can be expressed by magnitude time cos/sin angle
            se = sqrt(sx.^2 + sy.^2);
            phi = acos(sx ./ se);

            % The potential is given as a power series expansion
            %                     n   2m+3k
            % V    = sum    V    S   S       cos (3k*phi)
            %  nmk      nmk  nmk  a   e
            %
            % with n + 2m + 3k <= N. we use the expansion for N = 7.
            
            V = zeros(size(space.dvr{1,1}));
            
            for n = 0:7
                for m = 0:3
                    for k = 0:2
                        V = V + ...
                            obj.Vmat(n+1,m+1,k+1) * sa.^n .* se.^(2*m+3*k) .* cos(3*k*phi);
                    end
                end
            end
            
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            prt.error ('Code for calculation of forces still missing')
        end

        
    end
end