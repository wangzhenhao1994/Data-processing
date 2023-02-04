% Populations for each CHANNEL
function wave_pop (obj, psi, step)
global hamilt space

if ~strcmpi( obj.qua,'pop')
    prt.error ('Invalid choice of quantity for mean value / uncertainty')
end

for m = 1:hamilt.coupling.n_eqs
    obj.cha{m}(step) = sum ( abs(psi.dvr{m}(:)).^2 .* space.weight(:) );
end

end
