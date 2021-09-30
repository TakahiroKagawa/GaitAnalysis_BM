function[vel]=drv(pos,smpflq)
% drv: numerical derivative with central difference. 
% Input argument: time series of position data and sampling frequency.

dt = 1/smpflq; %time difference
[Datasize, Column ] = size(pos);
vel = zeros(Datasize,Column);
for cnt1 = 1:Column
    for cnt2 = 2:Datasize-1
        % central difference
        vel(cnt2,cnt1) = (pos(cnt2+1,cnt1)-pos(cnt2-1,cnt1))/(2*dt);
    end
    vel(1,cnt1) = vel(2,cnt1);
    vel(Datasize,cnt1) = vel(Datasize-1,cnt1);
end