
% main: main process of balance map analysis. 
% 1. File loading
% 2. Detect gait events
% 3. Preprocess for balance map
% 4. Paramters for balance map analysis
% 5. Boudnary of forward balance loss on balance map
% 6. Calculate states on the balance map
%% 1. File loading
clear;
fname='gait_sample.csv';

% import position data and contact data
data_table = readtable(fname);
data_array=table2array(data_table);

[DataSize, column] = size(data_array);
Time = data_array(:,1);
PosData_dmy = data_array(:,2:column-2);
LeftContact = data_array(:,column-1);
RightContact = data_array(:,column);
SamplingFreq = 1/(Time(2)-Time(1));
% create data structure of position data
PosDataLabel={'Hips';'Spine';'Spine1';'Neck';'Head';'LeftShoulder';'LeftArm';'LeftForeArm';'LeftHand';'RightShoulder';'RightArm';'RightForeArm';'RightHand';'LeftUpLeg';'LeftLeg';'LeftFoot';'LeftToe';'LeftHeel';'RightUpLeg';'RightLeg';'RightFoot';'RightToe';'RightHeel'};
SegmentSize = length(PosDataLabel);
for cnt1 = 1:SegmentSize
    PosData(cnt1).name = PosDataLabel{cnt1}; %data lable
    PosData(cnt1).data = PosData_dmy(:,3*cnt1-2:3*cnt1); % x y z data
    PosData(cnt1).vel = drv(PosData(cnt1).data,SamplingFreq); % x y z data
end
clear data_table data_array PosData_dmy;
%% 2. Detect gait events

%start from left heel strike
StepCnt=0;
State=1;
for cnt1 = 1:DataSize-1
    if State == 1 && LeftContact(cnt1)==0 && LeftContact(cnt1+1)==1
        StepCnt=StepCnt+1;
        LHS(StepCnt) = cnt1;
        State = 2;
    end
    if State == 2 && RightContact(cnt1)==1 && RightContact(cnt1+1)==0
        RTO(StepCnt) = cnt1;
        State = 3;
    end
    if State == 3 && RightContact(cnt1)==0 && RightContact(cnt1+1)==1
        RHS(StepCnt) = cnt1;
        State = 4;
    end
    if State == 4 && LeftContact(cnt1)==1 && LeftContact(cnt1+1)==0
        LTO(StepCnt) = cnt1;
        State = 1;
    end    
end

%% 3. Preprocess for balance map
param.Height=162;
param.Weight=52;
% Lowerleg, Upperleg, foot, forearm, upperarm, hand, head, Thorax, Abdomain, pelvis
param.SegmentMass=[0.0465, 0.1, 0.0145, 0.016, 0.028, 0.006, 0.081, 0.216, 0.139, 0.142]; 
% Lowerleg, Upperleg, foot, forearm, upperarm, hand, head, Thorax, Abdomain, pelvis
param.SegmentCOMLength=[0.567, 0.567, 0, 0.57, 0.564, 0, 1, 0.57, 0.564, 0.288];

% Calculate segment and body center of mass
[COMdata]=SegmentCOM(PosData,param); % calculate center of mass position of each body segment
COMSize=length(COMdata);
for cnt1 = 1:COMSize
    COMdata(cnt1).vel = drv(COMdata(cnt1).data,SamplingFreq);
end

% animation of stick picture with positions of segment and body center of mass 
% animation_COM(PosData, COMdata);

% variable to specify the parameter L1 and L2 from measured data
Lst=[];
Lsw=[];
% Segmentation of swing phase data
for cnt1 = 1:StepCnt-1
    
    % Right swing phase
    StartIdx = RTO(cnt1);
    EndIdx =RHS(cnt1)-1;
    StanceAnkle.data = PosData(16).data(StartIdx:EndIdx,2:3);
    StanceAnkle.vel = PosData(16).vel(StartIdx:EndIdx,2:3);
    Hip.data = PosData(1).data(StartIdx:EndIdx,2:3);
    Hip.vel = PosData(1).vel(StartIdx:EndIdx,2:3);
    x_st.data = (sum(param.SegmentMass(1:3))*COMdata(18).data(StartIdx:EndIdx,2:3)+(sum(2*param.SegmentMass(4:6))+sum(param.SegmentMass(7:10)))*COMdata(17).data(StartIdx:EndIdx,2:3))/(sum(param.SegmentMass(1:3))+2*sum(param.SegmentMass(4:6))+sum(param.SegmentMass(7:10)));
    x_st.data = x_st.data - StanceAnkle.data; 
    x_sw.data = COMdata(19).data(StartIdx:EndIdx,2:3) - Hip.data; 
    x_st.vel = (sum(param.SegmentMass(1:3))*COMdata(18).vel(StartIdx:EndIdx,2:3)+(sum(2*param.SegmentMass(4:6))+sum(param.SegmentMass(7:10)))*COMdata(17).vel(StartIdx:EndIdx,2:3))/(sum(param.SegmentMass(1:3))+2*sum(param.SegmentMass(4:6))+sum(param.SegmentMass(7:10)));
    x_st.vel = x_st.vel - StanceAnkle.vel; 
    x_sw.vel = COMdata(19).vel(StartIdx:EndIdx,2:3) - Hip.vel; 

    % Stance and Swing leg position and velocity
    STEP(2*cnt1-1).COMState = [x_st.data(:,1), x_sw.data(:,1), x_st.vel(:,1), x_sw.vel(:,1)];
    STEP(2*cnt1-1).Time = Time(StartIdx:EndIdx);
    
    for cnt2 = 1:length(x_st.data)
        Lst =[Lst;norm(x_st.data(cnt2,:))]; % L1
        Lsw =[Lsw;norm(x_sw.data(cnt2,:))]; % L2
    end
    
    
    % Left swing phase
    StartIdx = LTO(cnt1);
    EndIdx =LHS(cnt1+1);
    % Horizontal position/velocity of stance leg and swing leg
    StanceAnkle.data = PosData(21).data(StartIdx:EndIdx,2:3);
    StanceAnkle.vel = PosData(21).vel(StartIdx:EndIdx,2:3);
    Hip.data = PosData(1).data(StartIdx:EndIdx,2:3);
    Hip.vel = PosData(1).vel(StartIdx:EndIdx,2:3);
    x_st.data = (sum(param.SegmentMass(1:3))*COMdata(19).data(StartIdx:EndIdx,2:3)+(sum(2*param.SegmentMass(4:6))+sum(param.SegmentMass(7:10)))*COMdata(17).data(StartIdx:EndIdx,2:3))/(sum(param.SegmentMass(1:3))+2*sum(param.SegmentMass(4:6))+sum(param.SegmentMass(7:10)));
    x_st.data = x_st.data - StanceAnkle.data; 
    x_sw.data = COMdata(18).data(StartIdx:EndIdx,2:3) - Hip.data; 
    x_st.vel = (sum(param.SegmentMass(1:3))*COMdata(19).vel(StartIdx:EndIdx,2:3)+(sum(2*param.SegmentMass(4:6))+sum(param.SegmentMass(7:10)))*COMdata(17).vel(StartIdx:EndIdx,2:3))/(sum(param.SegmentMass(1:3))+2*sum(param.SegmentMass(4:6))+sum(param.SegmentMass(7:10)));
    x_st.vel = x_st.vel - StanceAnkle.vel; 
    x_sw.vel = COMdata(18).vel(StartIdx:EndIdx,2:3) - Hip.vel; 

    % Stance and Swing leg position and velocity
    STEP(2*cnt1).COMState = [x_st.data(:,1), x_sw.data(:,1), x_st.vel(:,1), x_sw.vel(:,1)];
    STEP(2*cnt1).Time = Time(StartIdx:EndIdx);

    for cnt2 = 1:length(x_st.data)
        Lst =[Lst;norm(x_st.data(cnt2,:))]; % L1
        Lsw =[Lsw;norm(x_sw.data(cnt2,:))]; % L2
    end
end


%% 4. Paramters for balance map analysis
% gravity acceleration
BMpar.g =9.8;

% sampling frequecy of mocap data
BMpar.SmpFrq = SamplingFreq;

% mass parameters
Weight = 52; % body weight
BMpar.MI = Weight * (sum(param.SegmentMass(1:3))+2*sum(param.SegmentMass(4:6))+sum(param.SegmentMass(7:10))); % inverted pendulum
BMpar.MS = Weight * sum(param.SegmentMass(1:3)); % simple pendulum

% link length: average of the length of measured data
BMpar.LI =mean(Lst); % inverted pendulum
BMpar.LS =mean(Lsw); % simple pendulum

% parameters of natural frequencies
BMpar.OmegaI = sqrt(BMpar.g/BMpar.LI); % inverted pendulum
BMpar.OmegaS = sqrt(BMpar.g/BMpar.LS); % simple pendulum

% paramter of mass ratio
coefK = BMpar.MI/(BMpar.MI+BMpar.MS);
% matrix of linear equation of motion
BMpar.A=[BMpar.OmegaI^2/coefK, (1-coefK)/coefK*BMpar.OmegaI^2;
    -BMpar.OmegaS^2/coefK, -BMpar.OmegaS^2/coefK];
[V,S]=eig(BMpar.A);

dummy = inv(V);
% matrix for diagonalization with normalization 
BMpar.T(1,1:2)=dummy(1,1:2)/(dummy(1,1)+dummy(1,2));
BMpar.T(2,1:2)=dummy(2,1:2)/(dummy(2,1)+dummy(2,2));
BMpar.Tinv=inv(BMpar.T);

% natural frequency after diagonaliztion
BMpar.ChiOmegaI = sqrt(S(1,1));
BMpar.ChiOmegaS = sqrt(-S(2,2));
% Analytical solution
% DeltaOmega = BMpar.OmegaI^2-BMpar.OmegaS^2
% BMpar.ChiOmegaI = sqrt(1/(2*coefK)*(DeltaOmega+sqrt(DeltaOmega^2+4*coefK*BMpar.OmegaI^2*BMpar.OmegaS^2)));
% BMpar.ChiOmegaS = sqrt(1/(2*coefK)*(-DeltaOmega+sqrt(DeltaOmega^2+4*coefK*BMpar.OmegaI^2*BMpar.OmegaS^2)));

% natural frequency ratio
BMpar.ChiOmega0 = BMpar.ChiOmegaI/BMpar.ChiOmegaS;

%% 5. Boudnary of forward balance loss on balance map

E0_Positive=5:-0.01:0.01;
E0_Negative=-0.001:-0.01:-0.75;

Delay_Positive=StabilityBoundaryPositive(E0_Positive,BMpar.ChiOmega0);
Delay_Negative=StabilityBoundaryNegative(E0_Negative,BMpar.ChiOmega0);

y_lim_p=length(E0_Positive);
yd_lim_p=length(E0_Positive);
for cnt1=1:length(E0_Positive)
    y_lim_p(cnt1)=sqrt(E0_Positive(cnt1))/BMpar.ChiOmega0*sinh(Delay_Positive(cnt1));
    yd_lim_p(cnt1)=sqrt(E0_Positive(cnt1))*cosh(Delay_Positive(cnt1));
end
y_lim_n=length(E0_Negative);
yd_lim_n=length(E0_Negative);
for cnt1=1:length(E0_Negative)
    y_lim_n(cnt1)=sqrt(-E0_Negative(cnt1))/BMpar.ChiOmega0*cosh(Delay_Negative(cnt1));
    yd_lim_n(cnt1)=sqrt(-E0_Negative(cnt1))*sinh(Delay_Negative(cnt1));
end
y_lim=[y_lim_p,y_lim_n];
yd_lim=[yd_lim_p,yd_lim_n];


%% 6. Calculate states on the balance map

for cnt1 = 1:length(STEP)
    % calculate the swing leg position from COM position of swing leg
    STEP(cnt1).COMState(:,2)=BMpar.LI/BMpar.LS*STEP(cnt1).COMState(:,2);
    STEP(cnt1).COMState(:,4)=BMpar.LI/BMpar.LS*STEP(cnt1).COMState(:,4);
end

for cnt1 = 1:length(STEP)

    [DataSize,~] = size(STEP(cnt1).COMState);
    % Calculate energy ratio, phase difference 
    [STEP(cnt1).EnergyRatio, STEP(cnt1).PhaseDifference, STEP(cnt1).ES, STEP(cnt1).Chi, STEP(cnt1).T] = EnergyPhase(STEP(cnt1).COMState, BMpar);
    
    % Calculate chi_0 and dot_chi_0
    for cnt2 = 1:length(STEP(cnt1).EnergyRatio)
        E0 = STEP(cnt1).EnergyRatio(cnt2);
        Phase0 = STEP(cnt1).PhaseDifference(cnt2);
        if E0>0
            STEP(cnt1).chi0(cnt2,1)=sqrt(E0)/BMpar.ChiOmega0*sinh(Phase0); %position
            STEP(cnt1).chi0(cnt2,2)=sqrt(E0)*cosh(Phase0); %velocity
        else
            STEP(cnt1).chi0(cnt2,1)=sqrt(-E0)/BMpar.ChiOmega0*cosh(Phase0); %position
            STEP(cnt1).chi0(cnt2,2)=sqrt(-E0)*sinh(Phase0); %velocity
        end
    end
    % calculate margin from boundaries of balance loss region
    STEP(cnt1).Margin = MarginFunction(STEP(cnt1).chi0,y_lim,yd_lim,BMpar);
end
%% Animation of stick picture and trajectory on balance map
% writerObj = VideoWriter([fname,'.avi']);
% open(writerObj);

ToeOff = [RTO(1),LTO(1),RTO(2),LTO(2),RTO(3),LTO(3),RTO(4),LTO(4),RTO(5),LTO(5)];
HeelStrike = [RHS(1),LHS(2),RHS(2),LHS(3),RHS(3),LHS(4),RHS(4),LHS(5),RHS(5),LHS(6)];
StartIdx = ToeOff(1);
EndIdx = HeelStrike(end);
for cnt1 =1 :length(PosData)
    PosData(cnt1).data = PosData(cnt1).data(StartIdx:EndIdx,:);
    PosData(cnt1).vel = PosData(cnt1).vel(StartIdx:EndIdx,:);
    PosMin(cnt1,:) = min(PosData(cnt1).data);
    PosMax(cnt1,:) = max(PosData(cnt1).data);

end

figure(1);clf;
BalanceMap(2,y_lim, yd_lim, BMpar.ChiOmega0);hold on;ax1=gca;
figure(1);set(gcf,'position', [50   500   960   420]);
s1=subplot(1,2,1);
figobj1=get(ax1,'children');
copyobj(figobj1,s1);axis([-0.6 0.6 0 1.2]);axis square;
axis([-1 1.5 -0.5 1.8])
xlabel('Position');
ylabel('Velocity');
set(gca,'Fontsize', 16);
close(2)

h_BalancemapPoint=animatedline('Marker','o','MarkerEdgeColor','b','MarkerSize',10);
h_Balancemapline=animatedline('LineStyle','-','Linewidth',2,'color','b');
subplot(1,2,2);
walk_bitsHeadTrunk=animatedline('LineStyle','-','Linewidth',2,'color','g');%Lineの設定
walk_bitsRarmbranch=animatedline('LineStyle','-','Linewidth',2,'color','b');
walk_bitsLarmbranch=animatedline('LineStyle','-','Linewidth',2,'color','r');
walk_bitsRfootbranch=animatedline('LineStyle','-','Linewidth',2,'color','b');
walk_bitsLfootbranch=animatedline('LineStyle','-','Linewidth',2,'color','r');



% definition of size of the graph
set(gca,'XLim',[min(PosMin(:,1)),max(PosMax(:,1))]);
set(gca,'YLim',[-0.5,1.0]);
set(gca,'ZLim',[0,1.8]);
axis equal;
grid on;
set(gca,'Fontsize', 16);
view(90,0)

DataSize = length(PosData(1).data);
ToeOff = ToeOff-(StartIdx-1)*ones(size(ToeOff));
HeelStrike = HeelStrike-(StartIdx-1)*ones(size(HeelStrike));
StepIdx = 1;
for cnt = 1:DataSize
   
    subplot(1,2,1);
    clearpoints(h_BalancemapPoint);
    if cnt>=ToeOff(StepIdx) && cnt<=HeelStrike(StepIdx)-1
        addpoints(h_BalancemapPoint,BMpar.ChiOmega0*STEP(StepIdx).chi0(cnt-ToeOff(StepIdx)+1,1),STEP(StepIdx).chi0(cnt-ToeOff(StepIdx)+1,2));
        addpoints(h_Balancemapline,BMpar.ChiOmega0*STEP(StepIdx).chi0(cnt-ToeOff(StepIdx)+1,1),STEP(StepIdx).chi0(cnt-ToeOff(StepIdx)+1,2));
    else
        if cnt > HeelStrike(StepIdx)
            StepIdx = StepIdx+1;
            clearpoints(h_Balancemapline);
        end
    end
    drawnow;

    subplot(1,2,2);
    clearpoints(walk_bitsHeadTrunk);
    for cnt1 = 1:5
        addpoints(walk_bitsHeadTrunk,PosData(cnt1).data(cnt,1),PosData(cnt1).data(cnt,2),PosData(cnt1).data(cnt,3));
    end
    
    clearpoints(walk_bitsLarmbranch); % clear line segment
    for cnt1 = 6:9
        addpoints(walk_bitsLarmbranch,PosData(cnt1).data(cnt,1),PosData(cnt1).data(cnt,2),PosData(cnt1).data(cnt,3));
    end
    
    clearpoints(walk_bitsRarmbranch);
    for cnt1 = 10:13
        addpoints(walk_bitsRarmbranch,PosData(cnt1).data(cnt,1),PosData(cnt1).data(cnt,2),PosData(cnt1).data(cnt,3));
    end
    
    clearpoints(walk_bitsLfootbranch);
    for cnt1 = [1, 14:18,16]
        addpoints(walk_bitsLfootbranch,PosData(cnt1).data(cnt,1),PosData(cnt1).data(cnt,2),PosData(cnt1).data(cnt,3));
    end
    
    clearpoints(walk_bitsRfootbranch);
    for cnt1 = [1,19:23,21]
        addpoints(walk_bitsRfootbranch,PosData(cnt1).data(cnt,1),PosData(cnt1).data(cnt,2),PosData(cnt1).data(cnt,3));
    end
    

    drawnow;
    % FrameMat=getframe(gcf);
    % writeVideo(writerObj,FrameMat);

end
%close(writerObj);
