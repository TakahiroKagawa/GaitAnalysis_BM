function [Position] = ToeHeelPosition( skeleton )
% Left Foot
Toe=[0, skeleton(20).offset(2),skeleton(20).offset(3)]';%ç∂í‹êÊ
Heel=[0, skeleton(25).offset(2),-0.05]';%ç∂í‹êÊ
FootAxis=[0, skeleton(25).offset(2),0]';%ë´éÒåvéZóp
LFootSkeletonID = 19;
RFootSkeletonID = 24;
DataSize = length(skeleton(1).Dxyz);
for cnt1 = 1:DataSize
    LeftToe(cnt1,:) = (skeleton(LFootSkeletonID).trans(:,:,cnt1)*[Toe;1])';
    LeftHeel(cnt1,:) = (skeleton(LFootSkeletonID).trans(:,:,cnt1)*[Heel;1])';
    RightToe(cnt1,:) = (skeleton(RFootSkeletonID).trans(:,:,cnt1)*[Toe;1])';
    RightHeel(cnt1,:) = (skeleton(RFootSkeletonID).trans(:,:,cnt1)*[Heel;1])';
end
Position = [-LeftToe(:,1),LeftToe(:,[3,2]),-LeftHeel(:,1),LeftHeel(:,[3,2]),-RightToe(:,1),RightToe(:,[3,2]),-RightHeel(:,1),RightHeel(:,[3,2])];

