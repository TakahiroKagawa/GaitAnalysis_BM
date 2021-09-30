function [COMdata]=SegmentCOM(PosData,param)
% SegmentCOM: center of mass positions of body segments 
% Input argument: position data structure and parameter structure.

%Spine-Hips
Prox=PosData(1).data; % Proximal Joint
Distal=PosData(2).data; % Distal joint
COMdata(1).name = 'Pelvis';
COMdata(1).data = Distal+(Prox-Distal)*param.SegmentCOMLength(10);
%Spine1-Spine
Prox=PosData(2).data; % Proximal Joint
Distal=PosData(3).data; % Distal joint
COMdata(2).name = 'Abdomain';
COMdata(2).data = Distal+(Prox-Distal)*param.SegmentCOMLength(9);
%Neck-Spine1
Prox=PosData(3).data; % Proximal Joint
Distal=PosData(4).data; % Distal joint
COMdata(3).name = 'Thorax';
COMdata(3).data = Distal+(Prox-Distal)*param.SegmentCOMLength(8);
%Head
COMdata(4).name ='Head';
COMdata(4).data = PosData(5).data;


%LeftArm-LeftForeArm
Prox=PosData(7).data; % Proximal Joint
Distal=PosData(8).data; % Distal joint
COMdata(5).name = 'LeftUpperarm';
COMdata(5).data = Distal+(Prox-Distal)*param.SegmentCOMLength(5);
%LeftForeArm-LeftHand
Prox=PosData(8).data; % Proximal Joint
Distal=PosData(9).data; % Distal joint
COMdata(6).name = 'LeftForearm';
COMdata(6).data = Distal+(Prox-Distal)*param.SegmentCOMLength(4);
%LeftHand
COMdata(7).name = 'LeftHand';
COMdata(7).data = PosData(9).data;
%RightArm-RightForeArm
Prox=PosData(11).data; % Proximal Joint
Distal=PosData(12).data; % Distal joint
COMdata(8).name = 'RightUpperarm';
COMdata(8).data = Distal+(Prox-Distal)*param.SegmentCOMLength(5);
%RightForeArm-RightHand
Prox=PosData(12).data; % Proximal Joint
Distal=PosData(13).data; % Distal joint
COMdata(9).name = 'RightForearm';
COMdata(9).data = Distal+(Prox-Distal)*param.SegmentCOMLength(4);
%RightHand
COMdata(10).name = 'RightHand';
COMdata(10).data = PosData(13).data;

%LeftUpLeg-LeftLeg
Prox=PosData(14).data; % Proximal Joint
Distal=PosData(15).data; % Distal joint
COMdata(11).name = 'LeftUpperleg';
COMdata(11).data = Distal+(Prox-Distal)*param.SegmentCOMLength(2);
%LeftLeg-LeftFoot
Prox=PosData(15).data; % Proximal Joint
Distal=PosData(16).data; % Distal joint
COMdata(12).name = 'LeftLowerleg';
COMdata(12).data = Distal+(Prox-Distal)*param.SegmentCOMLength(1);
%LeftFoot
COMdata(13).name = 'LeftFoot';
COMdata(13).data = PosData(16).data;
%RightUpLeg-RightLeg
Prox=PosData(19).data; % Proximal Joint
Distal=PosData(20).data; % Distal joint
COMdata(14).name = 'RightUpperleg';
COMdata(14).data = Distal+(Prox-Distal)*param.SegmentCOMLength(2);
%RightLeg-RightFoot
Prox=PosData(20).data; % Proximal Joint
Distal=PosData(21).data; % Distal joint
COMdata(15).name = 'RightLowerleg';
COMdata(15).data = Distal+(Prox-Distal)*param.SegmentCOMLength(1);
%RightFoot
COMdata(16).name = 'RightFoot';
COMdata(16).data = PosData(21).data;

% COM of the left leg
COMdata(18).name = 'Lleg_COM';
COMdata(18).data = (param.SegmentMass(1)*COMdata(12).data+param.SegmentMass(2)*COMdata(11).data+param.SegmentMass(3)*COMdata(13).data)/sum(param.SegmentMass(1:3));
% COM of the right leg
COMdata(19).name = 'Rleg_COM';
COMdata(19).data = (param.SegmentMass(1)*COMdata(15).data+param.SegmentMass(2)*COMdata(14).data+param.SegmentMass(3)*COMdata(16).data)/sum(param.SegmentMass(1:3));
% COM of body
R_arm_dmy = param.SegmentMass(4)*COMdata(9).data+param.SegmentMass(5)*COMdata(8).data+param.SegmentMass(6)*COMdata(10).data;
L_arm_dmy = param.SegmentMass(4)*COMdata(6).data+param.SegmentMass(5)*COMdata(5).data+param.SegmentMass(6)*COMdata(7).data;
Body_dmy = param.SegmentMass(7)*COMdata(4).data+param.SegmentMass(8)*COMdata(3).data+param.SegmentMass(9)*COMdata(2).data+param.SegmentMass(10)*COMdata(1).data;
COMdata(17).name = 'HeadArmTrunk_COM';
COMdata(17).data = (R_arm_dmy+L_arm_dmy+Body_dmy)/(2*sum(param.SegmentMass(4:6))+sum(param.SegmentMass(7:10)));
%COM of whole body
COMdata(20).name = 'Wholebody_COM';
COMdata(20).data = (sum(param.SegmentMass(1:3))*COMdata(18).data+sum(param.SegmentMass(1:3))*COMdata(19).data+(sum(param.SegmentMass(4:6))+sum(param.SegmentMass(4:10)))*COMdata(17).data)/(2*sum(param.SegmentMass(1:3))+2*sum(param.SegmentMass(4:6))+sum(param.SegmentMass(7:10)));