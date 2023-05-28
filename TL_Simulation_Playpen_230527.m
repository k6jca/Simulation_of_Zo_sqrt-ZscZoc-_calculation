% TL_Simulation_Playpen_230527
% Version 1.0
%
% k6jca
% 27 May 23
%
% Simulate calibrating a "virtual" vna with ideal or real standards,
% and then measure Zsc and Zoc of a simulated transmission line
% terminated with either ideal or "real" short and open.
%
% Then calculate Zo = sqrt(Zsc*Zoc) and plot.


clear all;
close all;
clc;

%**************************************************************************
%**************************************************************************
%
% ***   RUN OPTIONS   ***
%
% First we select what our "VNA" uses as the cal-standard characteristics.
% Note that these are used to calculated G1, G2, and G3 (the "ACTUAL"
% Gammas of the SOL standards).
% 
% 1 selects the ideal characteristics (so G1, G2, and G3 will be -1, +1, and 0)
% 0 selects the "real world" characteristics.

characterize_with_ideal_SO = 1;  % for both short and open
characterize_with_ideal_Ld = 1;  % for both short and open

% Next we select which version of the SOL we will "measure" for
% calibration.
%
% 1 selects the ideal standard (so Gm1, Gm2, and Gm3 will be -1, +1, and 0)
% 0 selects the "real world" standards.

cal_with_ideal_SO = 1;  % for both short and open
cal_with_ideal_Ld = 1;  % for both short and open

% 1 selects Coax terminations are "ideal" Short & Open.  
% 0 selects the "real world" Short and Open

terminate_line_with_ideal_SO = 0;

%
%
%**************************************************************************
%**************************************************************************


% define frequency start and stop
start_freq = 1;    % MHz
stop_freq = 100;   % MHz
d_freq = 0.1;      % MHz

num_pts = (stop_freq - start_freq)/d_freq + 1

Zo = 50; % ohms


% ***************************
%
% Define Test Conditions, e.g.
%   o  Cal Standards 
%   o  Transmission Line
%   o  Testing Frequency Range, System Zo
%
% ***************************
%

% **********
% Cal Standards Actual Characteristics
% **********

%  Short
c0 = 50e-15;     % Value used in original NanoVNA code.
% c0 = 1e-40;      % use in lieu of 0 so that there is no div by 0
c1 = 0;
c2 = 0;
c3 = 0;
open_offset_delay = 0;  %
open_offset_delay = 17.5e-12;  % See:https://k6jca.blogspot.com/2019/09/notes-on-compensation-of-vna.html

% Open
l0 = 0;
l1 = 0;
l2 = 0;
l3 = 0;
short_offset_delay = 0;
short_offset_delay = 17.8e-12;  % See: https://k6jca.blogspot.com/2019/09/notes-on-compensation-of-vna.html


% Load
% rld =               50;      %  0 dB Return Loss
% rld =                 50.3;  % 50 dB Return Loss
rld =               51;      % 40 dB Return Loss
% rld =               53;      % 31 dB Return Loss
% rld =               60;      % 21 dB Return Loss
ld_offset_delay =   0;
% ld_offset_delay =   100e-12;



% **********
% Transmission Line
% (length, Zo, alpha, beta)
%
% source: https://scikit-rf.readthedocs.io/en/latest/examples/networktheory/Transmission%20Line%20Properties%20and%20Manipulations.html
% **********

len = 5;         % in meters
Zo_tl = 75;
VF = 0.67;
alpha =  0.02;      % atten, Neper per meter

c = 299792458;    % speed of light, m/sec

% ****************************
% ****************************
%
% Step through Frequency Range
% and run frequency-dependent 
% calculations
%
% ****************************
% ****************************

for k = 1:num_pts
    freq(k) = (start_freq + (k-1)*d_freq)*1e6; % multiply by 1e6 to get MHz
    w(k) = 2*pi*freq(k);

    F = freq(k);

    % ***************************
    %
    % Calculate Gammas of Cal Standards
    %
    % ***************************
    %
    %  G1  = Actual Gamma of Standard 1 = SHORT
    %  G2  = Actual Gamma of Standard 2 = OPEN
    %  G3  = Actual Gamma of Standard 3 = LOAD
    %  Gm1 = Measured Gamma of Standard 1
    %  Gm2 = Measured Gamma of Standard 2
    %  Gm3 = Measured Gamma of Standard 3


    % First, calculate actual values of the standards
    Lshort(k) = l0 + l1*F + l2*F^2 + l3*F^3;
    Zshort(k) = 1i*w(k)*Lshort(k);
    Gshort(k) = z2gamma(Zshort(k),50);
    Gshort(k) = Gshort(k)*exp(-2*1i*(2*pi()*F*short_offset_delay)); % apply offset delay to Gamma
    Zshort(k) = gamma2z(Gshort(k),50);    % convert Gamma back into Z

    Copen(k) = c0 + c1*F + c2*F^2 + c3*F^3;
    Zopen(k) = -1i/(w(k)*Copen(k));
    Gopen(k) = z2gamma(Zopen(k),50);
    Gopen(k) = Gopen(k)*exp(-2*1i*(2*pi()*F*open_offset_delay)); % apply offset delay to Gamma
    Zopen(k) = gamma2z(Gopen(k),50);       % convert Gamma back into Z

    Zld(k) = rld; 
    Gld(k) = z2gamma(Zld(k),50);
    Gld(k) = Gld(k)*exp(-2*1i*(2*pi()*F*ld_offset_delay));  % apply offset delay to Gamma
    Zld(k) = gamma2z(Gld(k),50);      % convert Gamma back into Z


    % Second, define the actual Gammas of the standards being used for cal.
    % These will either be IDEAL or actual values.
    if characterize_with_ideal_SO  % for both short and open
        G1(k) = -1;
        G2(k) = 1;
        G3(k) = 0;
    else        
        G1(k) = (Zshort(k)-Zo)/(Zshort(k)+Zo);        
        G2(k) = (Zopen(k)-Zo)/(Zopen(k)+Zo);
        G3(k) = (Zld(k)-Zo)/(Zld(k)+Zo);
    end
    if characterize_with_ideal_Ld  % for both short and open
        G3(k) = 0;
    else        
        G3(k) = (Zld(k)-Zo)/(Zld(k)+Zo);
    end

    % Third, select which SOL is to be used for the "measured"
    % Gamma calculation.
    if cal_with_ideal_SO  % for both short and open
        Gm1(k) = -1    % short
        Gm2(k) = +1;    % open
        Gm3(k) = 0;
    else
        Gm1(k) = (Zshort(k)-Zo)/(Zshort(k)+Zo);
        Gm2(k) = (Zopen(k)-Zo)/(Zopen(k)+Zo);
        Gm3(k) = (Zld(k)-Zo)/(Zld(k)+Zo);
    end
    if cal_with_ideal_Ld  % for both short and open
        Gm3(k) = 0;
    else
        Gm3(k) = (Zld(k)-Zo)/(Zld(k)+Zo);
    end

     
    
    % ***************************
    %
    % Calculate S11 Error Correction Terms
    %
    % ***************************
    
    e00(k) = (G1(k)*G2(k)*Gm1(k)*Gm3(k) - G1(k)*G3(k)*Gm1(k)*Gm2(k) - G1(k)*G2(k)*Gm2(k)*Gm3(k) + G2(k)*G3(k)*Gm1(k)*Gm2(k) + ...
        G1(k)*G3(k)*Gm2(k)*Gm3(k) - G2(k)*G3(k)*Gm1(k)*Gm3(k))/(G1(k)*G2(k)*Gm1(k) - G1(k)*G2(k)*Gm2(k) - G1(k)*G3(k)*Gm1(k) +...
        G1(k)*G3(k)*Gm3(k) + G2(k)*G3(k)*Gm2(k) - G2(k)*G3(k)*Gm3(k));
     
     
    e11(k) = (G1(k)*Gm2(k) - G2(k)*Gm1(k) - G1(k)*Gm3(k) + G3(k)*Gm1(k) + G2(k)*Gm3(k) - G3(k)*Gm2(k))/...
        (G1(k)*G2(k)*Gm1(k) - G1(k)*G2(k)*Gm2(k) - G1(k)*G3(k)*Gm1(k) + G1(k)*G3(k)*Gm3(k) + G2(k)*G3(k)*Gm2(k)...
        - G2(k)*G3(k)*Gm3(k));
     
     
    d_e(k) = -(G1(k)*Gm1(k)*Gm2(k) - G1(k)*Gm1(k)*Gm3(k) - G2(k)*Gm1(k)*Gm2(k) + G2(k)*Gm2(k)*Gm3(k) + G3(k)*Gm1(k)*Gm3(k) - ...
        G3(k)*Gm2(k)*Gm3(k))/(G1(k)*G2(k)*Gm1(k) - G1(k)*G2(k)*Gm2(k) - G1(k)*G3(k)*Gm1(k) + G1(k)*G3(k)*Gm3(k) + ...
        G2(k)*G3(k)*Gm2(k) - G2(k)*G3(k)*Gm3(k));
     
  
    
    % ***************************
    %
    % Shorted Line:
    %   Given short (ideal or actual)
    %   o  Calculate "Measured Gamma"
    %   o  Calculate "Actual Gamma" from "Measured Gamma"
    %   o  Calculate Zin from Actual Gamma
    %
    % ***************************
    
    % Define transmission line beta and gamma
    beta(k) = (w(k)/c)*VF;
    gamma(k) = alpha + 1i*beta(k);

    % to calculate Zin, first need Gamma of the short terminating the
    % transmission line.

    if terminate_line_with_ideal_SO 
        Gamma_SHORTterm = (0-Zo_tl)/(0+Zo_tl);
    else
        Gamma_SHORTterm = (Zshort(k)-Zo_tl)/(Zshort(k)+Zo_tl);
    end
    % calculate Z looking into the transmission line
    zsc(k) = Zo_tl*(1+Gamma_SHORTterm*exp(-2*gamma(k)*len))/(1-Gamma_SHORTterm*exp(-2*gamma(k)*len));
    
    % Measure S11 of terminated line on "virtual vna".  Then error-correct
    % S11.  Finally convert the corrected S11 to corrected Zsc.
    S11_SHORTterm(k) = z2gamma(zsc(k),50);
    S11_SHORTterm_corrected(k) = (S11_SHORTterm(k)-e00(k))/(e11(k)*S11_SHORTterm(k)-d_e(k));
    zsc_corrected(k) = gamma2z(S11_SHORTterm_corrected(k),50);


    % ***************************
    %
    % Open Line:
    %   o  Calculate "Measured Gamma"
    %   o  Calculate "Actual Gamma" from "Measured Gamma"
    %   o  Calculate Zin from Actual Gamma
    %
    % ***************************
    
    % to calculate Zin, need first need Gamma of the open terminating the
    % transmission line.

    if terminate_line_with_ideal_SO 
        Gamma_OPENterm = (1e20 - Zo_tl)/(1e20 + Zo_tl);
    else
        Gamma_OPENterm = (Zopen(k)-Zo_tl)/(Zopen(k)+Zo_tl);
    end
    % calculate Z looking into the transmission line:
    zoc(k) = Zo_tl*(1+Gamma_OPENterm*exp(-2*gamma(k)*len))/(1-Gamma_OPENterm*exp(-2*gamma(k)*len));

    % Measure S11 of terminated line on "virtual vna".  Then error-correct
    % S11.  Finally convert the corrected S11 to a corrected Zoc.
    S11_OPENterm(k) = z2gamma(zoc(k),50);
    S11_OPENterm_corrected(k) = (S11_OPENterm(k)-e00(k))/(e11(k)*S11_OPENterm(k)-d_e(k));
    zoc_corrected(k) = gamma2z(S11_OPENterm_corrected(k),50);
    
    % ***************************
    %
    % Calculate Zo = SQRT(zsc*zoc)
    %
    % ***************************
        
    zo_calc(k) = sqrt(zsc_corrected(k)*zoc_corrected(k));


end


% ***************************
%
% Plot results
%
% ***************************

h=figure()
t = tiledlayout(2,1)
ax1 = nexttile
plot(ax1,freq/1e6,real(zo_calc),'b-'); 
title(ax1,'Real (Ro) of Zo')
grid(ax1,'on');
 grid(ax1,'minor');
% ylim(ax1, [0,150]);
% xlim(ax1, [0,100]);
% legend(ax1,'SimSmith, no error','Gamma(short)*exp(-j*0.01)','location','southwest');
% legend(ax1,'SimSmith, no error','SimSmith, reduced precision','location','northeast');
 

ax2 = nexttile
plot(ax2,freq/1e6,imag(zo_calc),'b-'); 
title(ax2,'Imaginary (Xo) of Zo')
grid(ax2,'on');
 grid(ax2,'minor');
% ylim(ax2, [-50,50]);
% xlim(ax2, [0,100]);
% legend(ax2,'SimSmith, no error','Gamma(short)*exp(-j*0.01)','location','southwest');
% legend(ax2,'SimSmith, no error','SimSmith, reduced precision','location','northeast');

if (characterize_with_ideal_SO)
    char_SO_text = ' Ideal';
else
    char_SO_text = ' Actual';
end
if (characterize_with_ideal_Ld)
    char_Ld_text = ' Ideal';
else
    char_Ld_text = ' Actual';
end
if (cal_with_ideal_SO)
    cal_SO_text = ' Ideal';
else
    cal_SO_text = ' Actual';
end
if (cal_with_ideal_Ld)
    cal_Ld_text = ' Ideal';
else
    cal_Ld_text = ' Actual';
end
if (terminate_line_with_ideal_SO)
    term_SO_text = ' Ideal';
else
    term_SO_text = ' Actual';
end
test_text1 = strcat('Short/Open Characteristics:',char_SO_text,...
    ', Short/Open Cal with:',cal_SO_text);
test_text2 = strcat('50 ohm Load Characteristics:',char_Ld_text,...
    ', 50 ohm Load Cal with:',cal_Ld_text);
test_text3 = strcat('Terminate Line with:', term_SO_text,' Short/Open');
title(t,{'Zo Calculation of Simulated 75 ohm coax where Zo = Ro+jXo = SQRT(Zoc*Zsc)',...
    'Simulation Info:',test_text1,test_text2, test_text3})
% title(t,'My Title')
xlabel(t,'Frequency (MHz)')
ylabel(t,'Ohms')

set(h, 'Position', [50,350,1050, 600]);


