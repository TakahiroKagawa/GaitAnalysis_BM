function [delta]=StabilityBoundaryPositive(E0_array,omega0)
%options=optimoptions(@fsolve,'Display','off');
delta=zeros(1,length(E0_array));

for cnt1=1:length(E0_array)
    %BinarySearch 
    deltaMax = 3;
    deltaMin = 0;
    T = [0:0.001:2];
    ipsi = 1;
    while ipsi>0.001
        %BinarySearch 
        e=sin(T)-sqrt(E0_array(cnt1))/omega0*sinh(omega0*T+(deltaMax+deltaMin)/2);
        if max(e)*min(e)<0 % existing zero crossing
            deltaMin=(deltaMax+deltaMin)/2;
        else
            deltaMax=(deltaMax+deltaMin)/2;
        end
        ipsi=deltaMax-deltaMin;
    end
    delta(cnt1)=deltaMin;
end


