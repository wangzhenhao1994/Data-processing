%------------------------------------------------------------------
% Obtain spectrum as Fourier transform of autocorrelation function
%------------------------------------------------------------------

function spectrum ( obj )
global hamilt


%% Symmetrize time series: S(-t) = S(t)*

% Length of original time series: not counting zero time
nn = length(obj.acf)-1;

% Negative and zero times: flipping up-down
acf2 = flipud(conj(obj.acf(1:nn+1)));

% Positive times (excluding last time step)
acf2 = [acf2; obj.acf(2:nn)];

%% Fourier transform: time -> energy

% Frequency grid
fmax = pi / obj.s_delta;
df = fmax / nn;
obj.freq = [-nn:nn-1]'*df;

%% Time -> Frequency/energy

% Fourier transform
obj.spec = fftshift ( ifft  ( acf2 ));

% Truncate frequencies
inside = find ( obj.freq>hamilt.truncate.e_min & obj.freq<hamilt.truncate.e_max);
obj.freq = obj.freq(inside);
obj.spec = obj.spec(inside);

% Normalize intensities
spec_max = max(abs(obj.spec) );
if spec_max>0
    obj.spec = obj.spec / spec_max;
end

prt.disp ( '***************************************************************')
prt.disp ( 'Spectrum as Fourier transform of autocorrelation')
prt.disp ( '***************************************************************')
prt.disp ( ' ')
prt.disp (['Length of (mirrored!) time series       : ' int2str(2*nn)])
prt.disp (['Time resolution of autocorrelation      : ' num2str(obj.s_delta)])
prt.disp (['Maximum frequency/energy of spectrum    : ' num2str(fmax)])
prt.disp (['Frequency/energy resolution of spectrum : ' num2str(df)])
prt.disp ( ' ')

end
