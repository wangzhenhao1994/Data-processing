% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008,2010 Ulf Lorenz
%
% see the README file for license details.

function log ( step )
global expect hamilt info space state time



%% Construct header from several text strings
if strcmpi(info.program,'qm_bound') % TISE: Loop over eigenstates
    info.header1 = ['n = ', ...
        int2str(step-1+hamilt.eigen.start)]; % Start counting  from zero (for bound states)
else    % qm_propa or hand-crafted scripts
    info.header1 = ['t = ', ...
        num2str(time.steps.m_grid(step), '%10.4f'), ...
        ', step = ', num2str(step-1)];
end

switch(lower(class(state)))
    
    case {'ket','rho'}
        
        info.header2 = ['E = ', num2str(real(state.energy), '%12.6f'), ...
            ', N = ', num2str(real(state.norm), '%10.4f')];

        prt.disp ('***************************************************************')
        prt.disp (info.header1);
        prt.disp (info.header2);
        prt.disp ('***************************************************************')
        prt.disp (' ')
        for ii = 1:size(state.y,2)
            prt.disp ([int2str(ii) '-th observable: ' num2str(state.y(step,ii))])
        end
        prt.disp (' ')
        
    otherwise
        
        
        if hamilt.coupling.n_eqs==1
            estr = 'E = ';
        else
            if strcmpi(hamilt.coupling.represent,'adi')
                estr = 'E_{adi} = ';
            elseif strcmpi(hamilt.coupling.represent,'dia')
                estr = 'E_{dia} = ';
            end
        end
        
        info.header2 = [estr, ...
            num2str(expect.total(step), '%12.6f'), ...
            ', N = ', ...
            num2str(expect.pop.tot(step), '%10.4f')];
        
        prt.disp ('***************************************************************')
        prt.disp (info.header1);
        prt.disp (info.header2);
        prt.disp ('***************************************************************')
        prt.disp (' ')
        
        %% Expectation values of individual channel wavefunctions
        for m=1:hamilt.coupling.n_eqs
            
            % If population exceeds threshold
            if expect.pop.cha{m}(step)>expect.min_pop
                
                % Population
                if hamilt.coupling.n_eqs>1
                    if strcmpi(hamilt.coupling.represent,'adi')
                        prt.disp (['Adiabatic state ',int2str(m),': Population ',num2str(expect.pop.cha{m}(step),'%10.4f')])
                    elseif strcmpi(hamilt.coupling.represent,'dia')
                        prt.disp (['Diabatic state ',int2str(m),' (',hamilt.coupling.labels{m},') : Population ',num2str(expect.pop.cha{m}(step),'%10.4f')])
                    end
                    prt.disp ('------------------------------------------')
                    prt.disp (' ')
                end
                
                % additional multiplicative operators
                if isfield(hamilt, 'amo')
                    for p=1:length(hamilt.amo)
                        if ~isempty (hamilt.amo{p})
                            if ~isempty (hamilt.amo{p}.dvr)
                                prt.disp ([hamilt.amo{p}.label ': ',num2str(expect.amo{p}.cha{m}(step),'%12.6f')])
                            end
                        end
                    end
                    prt.disp (' ')
                end
                
                % Position/momentum
                for k = 1:space.n_dim
                    if space.n_dim > 1
                        prt.disp (['Degree of freedom : ' space.dof{k}.label])
                    end
                    prt.disp (['Position : ',num2str(expect.pos{k}.cha{m}(step),'%12.6f'),' +/- ',num2str(expect.pos{k}.unc{m}(step),'%12.6f')])
                    prt.disp (['Momentum : ',num2str(expect.mom{k}.cha{m}(step),'%12.6f'),' +/- ',num2str(expect.mom{k}.unc{m}(step),'%12.6f')])
                    prt.disp (['Uncertainty product : ',num2str(expect.pos{k}.unc{m}(step).*expect.mom{k}.unc{m}(step),'%12.6f')])
                    prt.disp (' ')
                end
                
                % Potential/kinetic/sum energy
                prt.disp (['Potential energy : ',num2str(expect.pot.cha{m}(step),'%15.8f'),' +/- ',num2str(expect.pot.unc{m}(step),'%12.6f')])
                prt.disp (['Kinetic energy   : ',num2str(expect.kin.cha{m}(step),'%15.8f'),' +/- ',num2str(expect.kin.unc{m}(step),'%12.6f')])
                prt.disp (['Sum of energies  : ',num2str(expect.pot.cha{m}(step)+expect.kin.cha{m}(step),'%15.8f')])
                prt.disp (' ')
                
            end
        end
        
        %% Expectation values of total  wavefunctions
        if hamilt.coupling.n_eqs>1
            
            % Population
            prt.disp (['Total: Population ',num2str(expect.pop.tot(step),'%10.4f')])
            prt.disp ('------------------------------------------')
            prt.disp (' ')
            
            % Projection
            if isfield(hamilt, 'amo')
                for p=1:length(hamilt.amo)
                    if ~isempty (hamilt.amo{p})
                        if ~isempty (hamilt.amo{p})
                            prt.disp ([hamilt.amo{p}.label ' : ',num2str(expect.amo{p}.tot(step),'%12.6f')])
                        end
                    end
                end
                prt.disp (' ')
            end
            
            % Position/momentum
            for k = 1:space.n_dim
                if space.n_dim > 1
                    prt.disp (['Degree of freedom : ' space.dof{k}.label])
                end
                prt.disp (['Position : ',num2str(expect.pos{k}.tot(step),'%12.6f')])
                prt.disp (['Momentum : ',num2str(expect.mom{k}.tot(step),'%12.6f')])
                prt.disp (' ')
            end
            
            % Potential/kinetic/sum energy
            prt.disp (['Potential energy : ',num2str(expect.pot.tot(step),'%12.6f')])
            prt.disp (['Kinetic energy   : ',num2str(expect.kin.tot(step),'%12.6f')])
            prt.disp (['Sum of energies  : ',num2str(expect.pot.tot(step)+expect.kin.tot(step),'%12.6f')])
            prt.disp (' ')
            prt.disp (['Incl. couplings  : ',num2str(expect.total(step),'%12.6f')])
            prt.disp (' ')
            
        end
        
end