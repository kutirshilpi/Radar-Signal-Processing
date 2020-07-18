
% Ambiguity function of linear FM signal
% Ashikur Rahman
% 7/18/2020

clear
Beta=1000e6; % FM sweep
Tau=10e-6; % Pulse time
t=-Tau:Tau/100:Tau;
L=length(t);
d=40/(L*Tau);
Fd=-20/Tau:40/(L*Tau):20/Tau-40/(L*Tau); % Target Doppler shift

A=1-abs(t)/Tau;
B=Fd*Tau+Beta*t;

ix=1;
iy=1;
for a=-Tau:Tau/100:Tau

    for f=-20/Tau:40/(L*Tau):20/Tau-40/(L*Tau)
    A=1-abs(a)/Tau;
    B=f*Tau+Beta*a;

    Af(ix,iy)=A.*abs(sin(pi*B.*A)/(pi*A.*B));
    ix=ix+1;
    end
    ix=1;
    iy=iy+1;
end

% surf(Af)
surf(Fd*Tau,t/Tau,Af)




%---------------------------------------------%

