% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%
% see the README file for license details.

function retval = subsref(obj, index)

global space

if index.type == '.'
	switch index.subs
        % first public properties that are required from all grid classes
		case 'label'
			retval = obj.label;
		case 'dof'
			retval = obj.dof;
        case 'n_pts'
            retval = obj.l_max - abs(obj.m_0) + 1;
        case 'dvr_min'
            retval = -1;
        case 'dvr_max'
            retval = +1;
        case 'fbr_min'
            retval = abs(obj.m_0);
        case 'fbr_max'
            retval = obj.l_max;
        case 'kin_max'
            % energy is l(l+1)/2mR^2
            if obj.nokin
                retval = 0;
            else
                if ~isempty(obj.R_0)
                    Rmin = obj.R_0;
                else
                    Rmin = min(space.dof{obj.R_dof}.x_grid);
                end
                retval = obj.l_max * (obj.l_max+1) / (2*obj.mass*Rmin^2);
            end
        case 'nokin'
            retval = obj.nokin;

            
        % Next the Legendre-specific properties
		case 'mass'
			retval = obj.mass;
		case 'R_dof'
			retval = obj.R_dof;
		case 'R_0'
			retval = obj.R_0;
		case 'l_max'
			retval = obj.l_max;
		case 'm_0'
			retval = obj.m_0;
            
        % Next the weights and grids
        case 'weight'
            retval = obj.weight;
        case 'x_grid'
            retval = obj.x_grid;
        case 'p_grid'
            retval = obj.p_grid;

        % Finally the private properties
		case 'kin'
			retval = obj.kin;
		case 'kin_expo'
			retval = obj.kin_expo;
        case 'trafo2fbr'
            retval = obj.trafo2fbr;
        case 'trafo2dvr'
            retval = obj.trafo2dvr;
		otherwise
			prt.error( ['unknown field ' index.subs ' of dof.legendre requested']);
	end
else
	prt.error('Tried unsupported access form for class dof.legendre');
end
