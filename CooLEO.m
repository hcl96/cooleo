%clear;clc;
%load('rpdata.mat'); %loads RPA matrix containing columns of station, radius, and convective heat transfer coefficient
station_rpa=rpdata(:,1);contour_rpa=rpdata(:,2);hc_rpa=rpdata(:,3); %creates data columns
% define Rp1 properties
Cp=2000; %J/kgK

% define engine properties
Taw=3300; %K adiabatic wall temperature obtained by Tc*Recovery factor=3500*0.93
Trg=310; %K set regen bulk temperature 40 deg c / 310K
Tcin(1)=380; %coolant injected at 400K
mdot_film=0.5; %kg/s
%regen eta by parts
hr_chamber=210; %desired regen cooling efficiency
hr_nozzle=300;

k=15; %W/mK Inconel thermal conductivity 
t=0.5*0.001; %0.5 mm wall thickness
hrc=(k*hr_chamber)/(t*hr_chamber+k); %composite heat transfer coefficient from film coolant to regen coolant, accounts for wall loss
hrn=(k*hr_nozzle)/(t*hr_nozzle+k);
Q=0; %initiliaze Q, heat transfer
Tr(1)=Tcin(1);
Trcin(1)=310;

%coolant becomes supercritical at 660K & 315 PSI or at 100 mm

for s=1:size(hc_rpa)-1
    hg=hc_rpa(s)*1000;%reads rpa matrix and convert to W/m^2K
    Rc(s)=contour_rpa(s)*0.001;%converts radial station distance to meters
    if station_rpa(s) < 100
        alpha=hrc+hg; %sum of heat transfer coefficient
        beta=hrc*Trg+hg*Taw; %calculates beta coeff.,see manual for details
        L=sqrt((station_rpa(s+1)-station_rpa(s))^2+(contour_rpa(s+1)-contour_rpa(s))^2)*0.001;  %calculates differential station length
        Tcout(s,1)=(beta/alpha)+(Tcin(s)-(beta/alpha))*exp(1)^((-2*pi*Rc(s)*L*alpha)/(mdot_film*Cp)); %calculates differential temperature
        Tcin(s+1)=Tcout(s,1); %creates buffer for iterative addition
    else
        alpha=hrn+hg; %sum of heat transfer coefficient
        beta=hrn*Trg+hg*Taw*0.7;
        L=sqrt((station_rpa(s+1)-station_rpa(s))^2+(contour_rpa(s+1)-contour_rpa(s))^2)*0.001;  %converts station length to m
        Tcout(s,1)=(beta/alpha)+(Tcin(s)-(beta/alpha))*exp(1)^((-2*pi*Rc(s)*L*alpha)/(mdot_film*Cp));
        Tcin(s+1)=Tcout(s,1);
    end
    var_s=s;
    if Tcout(s) > 660
        break
    end
end

if max(Tcout)>660
    for s1=var_s:size(hc_rpa)
        Tcout(s1,1)=660;
        Q=Q+hg*(Taw-Tcout(s1,1))*pi*2*Rc(s)*L;
        if Q>250000*mdot_film
            Tcout(s1,1)=1000; %calculates if heat of vaporization is within limits
        end
    end
end

Tcout(var_s+1,1)=Tcout(var_s); %adjust matrix dimensions for plot

xlabel('Station (mm)');ylabel('Temperature (K)');
plot(station_rpa,Tcout)
hold on
plot(station_rpa,10*contour_rpa,'black')
plot(station_rpa,hc_rpa*100,'red')
xlabel('Station (mm)');ylabel('Film Coolant Temperature (K)');

%% regen cooling calculation

k=0.1;visc=1.05*10^-3;visc_w=4.04*10^-4;
mdot_total=mdot_film+.44;
Pr=(visc*Cp)/k;
K=(k/hr_nozzle)*0.0214*(Pr^0.4)*((visc/visc_w)^0.14)*((1/visc)^0.8)*((8*mdot_total)/pi)^0.8;
K=K^1.25;

Rt=min(contour_rpa)*0.001;%radius of throat

d=0;    %initialize 0 diameter
error=1000;     %initialize
while abs(error) > 10 %solve d by error bound less than 0.001mm
    error=1.2*(pi*(2*Rt+0.8*(d+0.0005))/(d+0.0005))-K*(d^-2.25);
    d=d+0.000001;
end

N=1.2*(pi*(2*Rt+0.8*(d+0.0005))/(d+0.0005));
d_mm=d*1000; %converts to mm

%delta P= f*(L/D)*(pv^2/2)
rho=805;%kg/m^3
V_regen=(mdot_total/rho)/(N*(pi*d^2)/4);
Re=(rho*V_regen*d)/visc;
L=station_rpa(64)*0.001; %converts to meters
e=(44*10^-6)/d;
f=0.04;
del_P=f*((L*0.5)/d)*(rho*V_regen^2)/2;
del_P_psi=0.000145038*del_P;

fprintf('%2.f%% Film Cooling \n', 100*(mdot_film/1.5));
fprintf('Max temperature in film coolant: %3.f K \n', max(Tcout));
fprintf('Regen diameter: %.2f mm \n', d_mm);
fprintf('Regen number of channels: %2.f \n', N);
fprintf('Regen Vmax: %3.f m/s \n', V_regen);
fprintf('Regen Re: %3.f \n', Re);
fprintf('Regen Delta P: %3.f psi \n', del_P_psi);
