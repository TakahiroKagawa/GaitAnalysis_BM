function [Margin]=MarginFunction(chi0,y_lim,yd_lim,BMpar)
% MarginFunction: Calcurating margin values from the boundaries of forward and backward balance loss regions. 
% Input argument: state on the balance map, position and velocity of the boundary of forward balance loss, and paratemers for balance map.

[ChiLength,~]=size(chi0);
% chi0 is regarded as multiplied value of natural frequency ratio
chi0(:,1) = BMpar.ChiOmega0*chi0(:,1);

[BorderLength]=length(y_lim);
%Initialization of margin with sufficiently large values
Margin(:,1)=100000*ones(ChiLength,1);
Margin(:,2)=100000*ones(ChiLength,1);
NearestNeighbor=zeros(ChiLength,2);
Distance = zeros(BorderLength,1);

for cnt1 = 1:BorderLength
    % the distance of the boundary line from origin, which is used to
    % determine whether the state is inside or outside from the balance loss region 
    NormLim(cnt1)=norm([BMpar.ChiOmega0*y_lim(cnt1),yd_lim(cnt1)]);
end


for cnt1 = 1:ChiLength
    if chi0(cnt1,2)>1
        % When dot_chi0 > 1, the minimum distance is the horizontal
        % distance
        Margin(cnt1,1)=-chi0(cnt1,1);
    else
        % the distance of the boundary line from origin, which is used to
        % determine whether the state is inside or outside from the balance loss region 
        normChi=norm(chi0(cnt1,:));
        
        % positive border
        % find nearest neighbor point on the boundary of forward balance loss
        for cnt2 = 1:BorderLength
            Distance(cnt2)=norm(chi0(cnt1,:)-[BMpar.ChiOmega0*y_lim(cnt2),yd_lim(cnt2)]);
        end
        [MinDist,Idx]=min(Distance);
        NearestNeighbor(cnt1,:)=[BMpar.ChiOmega0*y_lim(Idx),yd_lim(Idx)];
        % calculate distance
        % If state is in the balane loss region, the distance is negative
        if norm(chi0(cnt1,:))<norm(NearestNeighbor(cnt1,:));
            Margin(cnt1,1)=MinDist;
        else
            if chi0(cnt1,1) < 0
                Margin(cnt1,1)=MinDist;
            else
                Margin(cnt1,1)=-MinDist;
            end
        end
    end
    % negative border
    EI=0.5*chi0(cnt1,2)^2-0.5*chi0(cnt1,1)^2;
    if EI<0 && chi0(cnt1,1)<0
        Margin(cnt1,2)=-abs(chi0(cnt1,1)+chi0(cnt1,2))/sqrt(2);
    else
        Margin(cnt1,2)=abs(chi0(cnt1,1)+chi0(cnt1,2))/sqrt(2);
    end
end
cnt1;