function [split,hand1,hand2] = schroedinger1D(varargin)
%------------------------------------------------------------
% Call: [E,E_fun] = schroedinger1D(nstates,name,N,show,m1,m2)
%------------------------------------------------------------
% Input: nstates        --- Number of states
% Input: name           --- Potential file name
%                           Format: N lines, 2 colums,
%                         1st column x values in A,
%                         2nd column energies in K
% Input: N              --- Number of steps in x
% Input: show           --- plot option: 1/0 for yes/no
% Input: m1,m2          --- dimer atomic masses in amu
%------------------------------------------------------------
% Output: Energy Eigenvalues E and Eigenfunctions E_fun
%------------------------------------------------------------
% Discription: This function calculates the Eigen-
%              energies and Eigenfunctions for a given
%              one-dimensional potential by finite
%              difference approach.
%-----------------------------------------------------------
% Author: Andreas Hauser
%-----------------------------------------------------------
% Date:  17/5/2010 mod. 28/9/2011 for diatomics
%                  mod. 30/7/2012 to read energies in K
%                  mod. 03/8/2012 ro read energies in kJ/mol
%                  mod. 29.01.2023 to read energies in Hartree
%-----------------------------------------------------------
nstates = varargin{1};
name = varargin{2};
N = varargin{3};
show = logical(varargin{4});
m1 = varargin{5};
m2 = varargin{6};

% Define constants
hbar = (6.6260755e-34)/(2*pi);
amu = 1.66053886e-27;
h2j = 4.359748e-18 ;% Hartree to Joule
j2meV = 1/1.6e-22; %Joule to milli EV
j2cm = 1/h2j*219474.628;
cm2K = 1.4387863;
K2cm = 1/cm2K;
Navo = 6.02214129e23;

% Calculate reduced mass
m = m1*m2/(m1+m2)*amu;

% Read potential and interpolate
dummy = load(name);

%xd = (dummy(:,1))*1e-10;   % von A in m
xd = (dummy(:,1))*5.29177210903*1e-11; % from bohr to m
yd = (dummy(:,2))*h2j;
yd = yd - min(yd); % Potentialminium auf Null setzen
if yd(1)>0
  xmin = min(xd);%-1.5e-10;
  xmax = max(xd);

  x = linspace(xmin,xmax,N);
  y = spline(xd,yd,x);
  dx = x(2)-x(1);

  % Generate tridiagonal matrix by discretization
  V = diag(y);
  T = (diag(-ones(1,N)*2) + diag(ones(1,N-1),1) + diag(ones(1,N-1),-1))/dx^2;
  H = sparse(2*m*V/hbar^2 - T);

  [E_fun,E] = eigs(H,nstates,'SA');

  E = diag(E*hbar^2/(2*m)); % in J

  E_cm = E*j2cm;
  E_meV = E*j2meV;


  % Plots if desired
  if show == 1
      %title('Eigenstates for given 1D-potential','FontSize',14,'interpreter','latex');
      xlabel('z (m)','FontSize',14,'FontName','Times')
      set(gca,'XMinorTick','on');
      ylabel('Energy (cm-1)','FontSize',14,'FontName','Times')
      hold on;

      u = length(y);
      asymp = y(u)*j2cm;
      %sfactor = 5e-21; % Scaling factor for density (energy in J)
      %sfactor = 10; % Scaling factor for density (energy in meV)
      sfactor = 750; % Scaling factor for density (energy in cm-1)
      for k=1:nstates
          plot(x,sfactor*(E_fun(:,k)).^2+E_cm(k)-asymp,'b');
          hand2 = plot(x,ones(1,N)*E_cm(k)-asymp,'b','LineWidth',1.0);
      end
      hand1 = plot(x,y*j2cm-asymp,'k','LineWidth',1.2);
      ymin = min(y*j2cm) - abs(E_cm(1))/2;
      ymax = max(E_cm) + abs(E_cm(1));
      axis([xmin xmax ymin-asymp ymax-asymp+100]);
  end

  %E_cm %-asymp
  split = []
  for k = 2:1:nstates
      split(end+1) = E_cm(k-1)-E_cm(k);
  end
else
  hand1 = 0;
  hand2 = 0;
  split = 0;
end
