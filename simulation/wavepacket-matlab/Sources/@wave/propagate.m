%--------------------------------------------------------------------------
%
% Propagate wavefunction(s) subject to a given Hamiltonian
% using one of the user-specified PDE solvers for TDSEs
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008,2010 Ulf Lorenz
%
% see the README file for license details.

function propagate ( obj, step )
global hamilt time

% Use one of the propagator classes of your choice; initialize ACF
if step==1
    time.steps.offset = + 1;
    init ( time.propa, obj );
    prt.disp ('***************************************************************')
    prt.disp ('Numeric propagation scheme:                       ')
    disp ( time.propa )
    prt.disp (' ')
    time.steps.acf    = zeros(length(time.steps.s_grid),1);
    time.steps.acf(1) = 1; % should hold by definition of ACF
    
else
    time.steps.offset = (step-2) * time.steps.s_number + 1;
    class_name = lower(class(time.propa)); 
    switch class_name (10:end)
        case {'diff_2','trotter','strang'}  % propagators that permit sub-stepping
            
            % Loop over substeps: Propagate and get autocorrelation
            for k = 1:time.steps.s_number
                propa ( time.propa, obj, k );
                time.steps.acf(k+time.steps.offset) = wave.braket(obj.ini,obj.dvr);
            end
            
        case {'cheby_real','cheby_imag'} % propagators that do not permit sub-stepping
            propa ( time.propa, obj );
    end
    
    % Apply negative imaginary potential (absorbing boundary conditions)
    % using a Trotter splitting
    if isfield (hamilt,'nip')
        for m=1:hamilt.coupling.n_eqs
            if ~isempty (hamilt.nip{m}.dvr)
                obj.dvr{m} = obj.dvr{m} .* exp( - hamilt.nip{m}.dvr ) ;
            end
        end
    end
    
end




