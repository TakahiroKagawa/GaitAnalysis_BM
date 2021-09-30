function [E0, Phase0, ES, chi, T] = EnergyPhase(COMstate, BMpar)
% EnergyPhase: Calcurating Energy ratio, phase difference from stance leg position and swing leg position and paramters for balancemap.
% In addition this function output energy of swing leg, non-dimensional position and time.
% Input argument: Swing leg position and Stance leg position and Paramters for balance map

DataSize=length(COMstate);
% variable transformation for diagonalization
x_prime(:,1)=BMpar.T(1,1)*COMstate(:,1)+BMpar.T(1,2)*COMstate(:,2);
x_prime(:,2)=BMpar.T(2,1)*COMstate(:,1)+BMpar.T(2,2)*COMstate(:,2);
x_prime(:,3)=BMpar.T(1,1)*COMstate(:,3)+BMpar.T(1,2)*COMstate(:,4);
x_prime(:,4)=BMpar.T(2,1)*COMstate(:,3)+BMpar.T(2,2)*COMstate(:,4);

for cnt1 = 1:DataSize
    % energy of stance and swing leg
    EI(cnt1) = (x_prime(cnt1,3)^2 - (BMpar.ChiOmegaI*x_prime(cnt1,1))^2)/2;
    ES(cnt1) = (x_prime(cnt1,4)^2 + (BMpar.ChiOmegaS*x_prime(cnt1,2))^2)/2;
    
    % phase of stance and swing leg
    if EI(cnt1)>0
        PhaseI(cnt1)=atanh(BMpar.ChiOmegaI*x_prime(cnt1,1)/x_prime(cnt1,3));
    else
        PhaseI(cnt1)=atanh(x_prime(cnt1,3)/(BMpar.ChiOmegaI*x_prime(cnt1,1)));
    end
    PhaseS(cnt1)=atan2(BMpar.ChiOmegaS*x_prime(cnt1,2), x_prime(cnt1,4));
    
    T(cnt1)=PhaseS(cnt1); % non-dimensional time
    E0(cnt1)=EI(cnt1)/ES(cnt1); % energy ratio
    Phase0(cnt1)=PhaseI(cnt1)-BMpar.ChiOmega0*PhaseS(cnt1); % phase difference
    
    % non-dimensional position
    if E0(cnt1)>0
        chi(cnt1,1)=sqrt(E0(cnt1))/BMpar.ChiOmega0*sinh(BMpar.ChiOmega0*T(cnt1)+Phase0(cnt1));
        chi(cnt1,3)=sqrt(E0(cnt1))*cosh(BMpar.ChiOmega0*T(cnt1)+Phase0(cnt1));
    else
        if x_prime(cnt1,1)>0
            chi(cnt1,1)=sqrt(-E0(cnt1))/BMpar.ChiOmega0*cosh(BMpar.ChiOmega0*T(cnt1)+Phase0(cnt1));
            chi(cnt1,3)=sqrt(-E0(cnt1))*sinh(BMpar.ChiOmega0*T(cnt1)+Phase0(cnt1));
        else
            chi(cnt1,1)=-sqrt(-E0(cnt1))/BMpar.ChiOmega0*cosh(BMpar.ChiOmega0*T(cnt1)+Phase0(cnt1));
            chi(cnt1,3)=-sqrt(-E0(cnt1))*sinh(BMpar.ChiOmega0*T(cnt1)+Phase0(cnt1));
        end
    end
    chi(cnt1,2)=sin(T(cnt1));
    chi(cnt1,4)=cos(T(cnt1));  
    
%     if E0(cnt1)>0
%         chi(cnt1,1)=sqrt(E0(cnt1))/BMpar.ChiOmega0*sinh(Phase0(cnt1));
%         chi(cnt1,3)=sqrt(E0(cnt1))*cosh(Phase0(cnt1));
%     else
%         if x_prime(cnt1,1)>0
%             chi(cnt1,1)=sqrt(-E0(cnt1))/BMpar.ChiOmega0*cosh(Phase0(cnt1));
%             chi(cnt1,3)=sqrt(-E0(cnt1))*sinh(Phase0(cnt1));
%         else
%             chi(cnt1,1)=-sqrt(-E0(cnt1))/BMpar.ChiOmega0*cosh(Phase0(cnt1));
%             chi(cnt1,3)=-sqrt(-E0(cnt1))*sinh(Phase0(cnt1));
%         end
%     end
end
%  [chi_h, ES] = x2chi(COMstate,BMpar);

end


function x_prime = x2prime(x,par)
% coefK = par.MI/(par.MI+par.MS);
% Sigma1 = par.ChiOmegaI^2;
% Sigma2 = -par.ChiOmegaS^2;
% % Alpha = 1/sqrt(par.OmegaI^4+(par.OmegaS^2+coefK*Sigma2)^2);
% % Beta = 1/sqrt(par.OmegaI^4+(par.OmegaS^2+coefK*Sigma1)^2);
% Alpha = 1/(par.OmegaI^2+(par.OmegaS^2+coefK*Sigma2)*par.LS/par.LI);
% Beta = 1/(par.OmegaI^2+(par.OmegaS^2+coefK*Sigma1)*par.LS/par.LI);
% P=[Alpha*par.OmegaI^2, Alpha*(par.OmegaS^2+coefK*Sigma2)*par.LS/par.LI;
% Beta*par.OmegaI^2, Beta*(par.OmegaS^2+coefK*Sigma1)*par.LS/par.LI];

x_prime(:,1)=par.T(1,1)*x(:,1)+par.T(1,2)*x(:,2);
x_prime(:,2)=par.T(2,1)*x(:,1)+par.T(2,2)*x(:,2);
x_prime(:,3)=par.T(1,1)*x(:,3)+par.T(1,2)*x(:,4);
x_prime(:,4)=par.T(2,1)*x(:,3)+par.T(2,2)*x(:,4);
end

function [chi, ES] = x2chi(x,par)
x_prime = x2prime(x,par);
[DataSize,~]=size(x_prime);
ES = zeros(DataSize,1);
chi = zeros(DataSize,4);
for cnt1 = 1:DataSize
    ES(cnt1) = 1/2*(x_prime(cnt1,4)^2)+1/2*(par.ChiOmegaS*x_prime(cnt1,2))^2;
    
    % State variables to Normalized position
    chi(cnt1,1)=par.ChiOmegaS/sqrt(2*ES(cnt1))*x_prime(cnt1,1);
    chi(cnt1,2)=par.ChiOmegaS/sqrt(2*ES(cnt1))*x_prime(cnt1,2);
    chi(cnt1,3)=1/sqrt(2*ES(cnt1))*x_prime(cnt1,3);
    chi(cnt1,4)=1/sqrt(2*ES(cnt1))*x_prime(cnt1,4);
end
end
