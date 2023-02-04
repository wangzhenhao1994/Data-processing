% Copyright (C) 2009-2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.


% Calculate the pump-probe signal for two-photon interferometry.
% We vary the timesteps between pump and probe in intervals of 10 fs,
% For each delay time, we run a calculation that loads the wavefunction
% from the pump-only simulation, and propagates it for 200 fs. Then, we
% calculate the observable. In our case, this is just the population of
% the 2^1\Pi_g state after the probe.
%
% Note that the reference (J.Cem.Phys. 100:5448) discusses in length the
% population caused by different processes (two photons taken from the
% pump pulse, probe pulse, or one from each). To somewhat disentangle the
% contributions to the final population, we propagate the pump result for
% the Sigma_g, Sigma_u, and Pi_g state separately. For non-overlapping
% pulses, this gives approximately the three wavefunctions psi.{11},
% psi.{22}, and psi.{21}. We can then overlap them as we want, and obtain
% the single contributions. See the reference for more details.
%
% The time delay is passed to the qm_init.m script via the global variable
% delay, the initial coefficients via init.coeffs.

global state

tp1    = zeros(201,1);     % ~2-photon absorbtion from pump pulse
tp2    = zeros(201,1);     % ~2-photon absorbtion from probe pulse
pp     = zeros(201,1);     % pump-probe signal
tpi    = zeros(201,1);     % two-photon interference
pptpi  = zeros(201,1);     % the pump-probe two-photon interference
tpppi  = zeros(201,1);     % the two-photon pump-probe interference
ii = 1;

% for practical reasons, we use the delay time in units of 10 fs.
for delay = 0:1:200
    
    qm_setup();
    qm_init(delay, 1);
    qm_propa();
    d22 = state.dvr{3};
    qm_cleanup();

    qm_setup();
    qm_init(delay, 3);
    qm_propa();
    d11 = state.dvr{3};
    qm_cleanup();

    qm_setup();
    qm_init(delay, 2);
	qm_propa();
    d21 = state.dvr{3};

    % calculate all the single terms
    global space
    tp1(ii) = sum(abs(d11(:)).^2 .* space.weight(:));
    tp2(ii) = sum(abs(d22(:)).^2 .* space.weight(:));
    pp(ii)  = sum(abs(d21(:)).^2 .* space.weight(:));

    tpi(ii)   = 2 * sum( real(conj(d11(:)) .* d22(:)) .* space.weight(:));
    pptpi(ii) = 2 * sum( real(conj(d21(:)) .* d22(:)) .* space.weight(:));
    tpppi(ii) = 2 * sum( real(conj(d21(:)) .* d11(:)) .* space.weight(:));
    
    qm_cleanup();

    ii = ii + 1;
end

total = tp1 + tp2 + pp + tpi + pptpi + tpppi;

figure(71);
clf;
plot(0:10:2000, total, 'LineWidth', 2);
xlabel('delay time (fs)');
ylabel('total Population of 2^1\Pi_g');
saveas(gcf, 'total.jpg');

figure(72);
clf;
plot(0:10:2000, tp1 + tp2, 'LineWidth', 2);
xlabel('delay time (fs)');
ylabel('Population of 2^1\Pi_g (both photons from one pulse)');
saveas(gcf, 'two-photon.jpg');

figure(73);
clf;
plot(0:10:2000, pp, 'LineWidth', 2);
xlabel('delay time (fs)');
ylabel('Population of 2^1\Pi_g (Pump-Probe signal)');
saveas(gcf, 'pp.jpg');

figure(74);
clf;
plot(0:10:2000, tpi, 'LineWidth', 2);
xlabel('delay time (fs)');
ylabel('Population of 2^1\Pi_g (two-photon interference)');
saveas(gcf, 'tpi.jpg');

figure(75);
clf;
plot(0:10:2000, pptpi, 'LineWidth', 2);
xlabel('delay time (fs)');
ylabel('Population of 2^1\Pi_g (interference between pump-probe and probe-probe)');
saveas(gcf, 'pptpi.jpg');

figure(76);
clf;
plot(0:10:2000, tpppi, 'LineWidth', 2);
xlabel('delay time (fs)');
ylabel('Population of 2^1\Pi_g (interference between pump-pump and pump-probe)');
saveas(gcf, 'tpppi.jpg');
