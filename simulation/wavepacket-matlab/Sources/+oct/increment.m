function delta = increment (iter)
%--------------------------------------------------------------------------
% Increment of total control functional
% see DOI:10.1063/1.1650297
%--------------------------------------------------------------------------
global control state

% Equation (17)
delta = control.optimal.alpha * trapz ...
    (control.t.steps, ...
    (2/control.optimal.zeta-1) * ...
    (control.u.backward - control.u.forward).^2 + ...
    (2/control.optimal.eta-1) * ...
    (control.u.backward - control.u.previous).^2 ./ ...
    control.u.shaping(1,:) );

% Equation (A5)
if ~isempty(state.D)
    dxT = control.x.previous(:,end) - control.x.forward(:,end);
    delta = delta + dot( dxT, state.D{control.optimal.terminal} * dxT);
end

% Case (C2), i.e. for D=|c><c|
if ~isempty(state.C) && state.Q{control.optimal.terminal}
    dxT = control.x.previous(:,end) - control.x.forward(:,end);
    delta = delta + abs ( state.C{control.optimal.terminal} * dxT )^2;
end

% Output
delta = math.real (delta);
prt.disp (['After ' int2str(iter) ' iterations:     DJ = ' num2str(delta)])

end


