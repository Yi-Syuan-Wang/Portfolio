%% Clear All
clc
close all
clear

bwI = 2;
skiI = 3;


%% Define the path & Import the Data
ID = 'T1';
Session = 'Session 1';
Task = 'P5'; 

codePath = 'D:\YS.W\Documents\NTU\Master Thesis\Code';
addpath(codePath);
dataPath = ['D:\YS.W\Documents\NTU\Master Thesis\Data','\',ID,'\',Session];
outputPath = 'D:\YS.W\Documents\NTU\Master Thesis\Data';
cd(dataPath);

load(['Rawdata_',Task,'.mat'])
Fx1 = Rawdata.ACHANNEL.ForceFx1;
Fy1 = Rawdata.ACHANNEL.ForceFy1;
Fz1 = Rawdata.ACHANNEL.ForceFz1;
Mx1 = Rawdata.ACHANNEL.MomentMx1;
My1 = Rawdata.ACHANNEL.MomentMy1;
Mz1 = Rawdata.ACHANNEL.MomentMz1;
Fx2 = Rawdata.ACHANNEL.ForceFx2;
Fy2 = Rawdata.ACHANNEL.ForceFy2;
Fz2 = Rawdata.ACHANNEL.ForceFz2;
Mx2 = Rawdata.ACHANNEL.MomentMx2;
My2 = Rawdata.ACHANNEL.MomentMy2;
Mz2 = Rawdata.ACHANNEL.MomentMz2;

Com = Rawdata.XYZPOS.CentreOfMass(:,3);
%Joint = struct('Ang', {Rawdata.XYZPOS.RAnkleAngles(:,1) Rawdata.XYZPOS.RKneeAngles(:,1) Rawdata.XYZPOS.RHipAngles(:,1)}); % Joint.Ang = {Ankle, Knee, Hip}
%Segment = struct('Axis',{Rawdata.XYZPOS.RTOO cat(3,Rawdata.XYZPOS.RFOO, Rawdata.XYZPOS.RFOA, Rawdata.XYZPOS.RFOL, Rawdata.XYZPOS.RFOP)... % Segment.Axis = {Toe(O), Foot(O,A,L,P), Tib, Fem}
    %cat(3,Rawdata.XYZPOS.RTIO, Rawdata.XYZPOS.RTIA, Rawdata.XYZPOS.RTIL, Rawdata.XYZPOS.RTIP)...
    %cat(3,Rawdata.XYZPOS.RFEO, Rawdata.XYZPOS.RFEA, Rawdata.XYZPOS.RFEL, Rawdata.XYZPOS.RFEP)},...
    %'F',{[Rawdata.ACHANNEL.ForceFx2 Rawdata.ACHANNEL.ForceFy2 Rawdata.ACHANNEL.ForceFz2] [] [] []},...
    %'M',{[Rawdata.ACHANNEL.MomentMx2 Rawdata.ACHANNEL.MomentMy2 Rawdata.ACHANNEL.MomentMz2] [] [] []});
Ang = [Rawdata.XYZPOS.RHipAngles(:,1) Rawdata.XYZPOS.RKneeAngles(:,1) Rawdata.XYZPOS.RAnkleAngles(:,1)]; % Hip flexion, knee flexion and ankle dorsiflexion as positive
Fem = cat(3,Rawdata.XYZPOS.RFEO, Rawdata.XYZPOS.RFEA, Rawdata.XYZPOS.RFEL, Rawdata.XYZPOS.RFEP);
Tib = cat(3,Rawdata.XYZPOS.RTIO, Rawdata.XYZPOS.RTIA, Rawdata.XYZPOS.RTIL, Rawdata.XYZPOS.RTIP);
Foot = cat(3,Rawdata.XYZPOS.RFOO, Rawdata.XYZPOS.RFOA, Rawdata.XYZPOS.RFOL, Rawdata.XYZPOS.RFOP);

%% Define Fixed Variables
fc = 10;
fsVi = 120;
fsFp = 960;
fsChange = fsVi/fsFp;
skiVi = skiI*fsVi;
skiFp = skiI*fsFp;
unitChange = 1/1000;
gap = 500;
fraFp = length(Fx1);
fraVi = length(Com); % Matrix size = m frams * n axes
dim = 3;
%% Filter & Smooth, unitChange
F1 = [-Fx1 -Fy1 -Fz1];
M1 = [-Mx1 -My1 -Mz1];
F2 = [-Fx2 -Fy2 -Fz2];
M2 = [-Mx2 -My2 -Mz2];
tFz = -(Fz1 + Fz2);
tFz = FilSmoo(tFz,2,fc,fsFp);
F1 = FilSmoo(F1,2,fc,fsFp);
M1 = FilSmoo(M1,2,fc,fsFp)*unitChange;
F2 = FilSmoo(F2,2,fc,fsFp);
M2 = FilSmoo(M2,2,fc,fsFp)*unitChange; % unit: Nmm -> Nm
Com = FilSmoo(Com,2,fc,fsVi,'low')/10;
Ang = FilSmoo(Ang,2,fc,fsVi,'low');

Fem = Fem*unitChange;
Tib = Tib*unitChange;
Foot = Foot*unitChange;

%% Calculate Anthropometric Data
BW = mean(tFz(fsFp*bwI:fsFp*bwI+fsFp)); 
stdBW = std(tFz(fsFp*bwI:fsFp*bwI+fsFp));
BM = BW/9.8;
femMass = 0.1*BM;
tibMass = 0.0465*BM;
footMass = 0.0145*BM;
femL = sqrt(sum((Fem(skiVi,:,4)-Fem(skiVi,:,1)).^2));
tibL = sqrt(sum((Tib(skiVi,:,4)-Tib(skiVi,:,1)).^2));
footL = sqrt(sum((Foot(skiVi,:,4)-Foot(skiVi,:,1)).^2));
femI = [femMass*(0.267*femL)^2 femMass*(0.267*femL)^2 femMass*(0.121*femL)^2]; % I = mx^2
tibI = [tibMass*(0.275*tibL)^2 tibMass*(0.281*tibL)^2 tibMass*(0.114*tibL)^2]; 
footI = [footMass*(0.124*footL)^2 footMass*(0.257*footL)^2 footMass*(0.245*footL)^2];

%% Define Eccentric Phase (eccSta, eccStaVi)
eccSta = find(tFz(skiFp:end) < BW-5*stdBW,1);  %When Fz over meanBW+-5*std, move forward 30ms
eccSta = eccSta+skiFp-round(30*fsFp*unitChange)-1;
eccStaVi = round(eccSta*fsChange); %Change the frame from FP to Vicon

%% Calculate Acceleration, Velocity and Displacement of COM from FP
comAFp = (tFz -BW)/(BW/10); % Unit of acceleartion: kg⋅m/s^2/kg = m/s^2
comVFp = zeros(length(comAFp),1); % Velocity before onset should be 0
comVFp(eccSta:end) = cumtrapz(comAFp(eccSta:end))/fsFp; % Use the trapezoid rule to integrate acceleration, unit of velocity: m/s
[minComVFp, minComVFpFra] = min(comVFp(eccSta:eccSta+fsFp));
minComVFpFra = minComVFpFra+eccSta-1;
unSta = minComVFpFra;
comDFp = zeros(length(comVFp),1); % Displacement before onset should be 0
comDFp(eccSta:end)= cumtrapz(comVFp(eccSta:end))/fsFp; % Use the trapezoid rule to integrate acceleration, unit of COM from FP: m

%% Define Concentric Phase
conSta = find(comVFp(minComVFpFra:end)>0, 1);
conSta = conSta+minComVFpFra-1;
conStaVi = round(conSta*fsChange);

%% Define Flight Phase
[maxComFp, maxComFpFra] = max(comDFp);
stdFlyF = std(tFz(maxComFpFra-150:maxComFpFra+150));
conEnd = find(tFz(conSta:end) <= 5*stdFlyF,1)+conSta-1;
conEndVi = round(conEnd*fsChange);

%% Duration of Each Phase
eccT = (conSta-eccSta)/fsFp;
conT = (conEnd-conSta+1)/fsFp;

%% Calculate JH & Depth
restCom = mean(Com(fsVi:fsVi*2));
[maxCom, maxComFra] = max(Com(conEndVi:conEndVi+60));
JH = maxCom - restCom;
[minCom,~] = min(Com(eccStaVi:conEndVi));
Depth = restCom - minCom;
%% Find Max. Joint Angle, ROM, Joint Angular Velocity
maxAng= Ang(conStaVi-1,:);
restAng = Ang(eccStaVi,:);
ROM = maxAng-restAng;
V = vel_diff(Ang,1/fsVi);
A = acc_diff(Ang,1/fsVi);
maxV = [max(V(eccStaVi:conStaVi-1,:));min(V(conStaVi:conEndVi,:))];

%% Print Kinematic Data
fprintf('Jump hiehgt=%.2f(cm) Depth =%.2f(cm)\nEccentric phase=%.4f(s)\t Concentric phase=%.4f(s)\n',JH,Depth,eccT,conT)
fprintf('Max. angle (eccentric) : Hip=%.4f(°)\t Knee=%.2f(°)\t Ankle=%.4f(°)\n',maxAng)
fprintf('ROM (eccentric) : Hip=%.4f(°)\t Knee=%.4f(°)\t Ankle=%.4f(°)\n',ROM)
fprintf('Max. angular velocity (eccentric) : Hip=%.4f(°/s)\t Knee=%.4f (°/s)\t Ankle=%.4f(°/s)\n',maxV(1,:))
fprintf('Max. angular velocity (concentric) : Hip=%.4f(°/s)\t Knee=%.4f(°/s)\t Ankle=%.4f(°/s)\n',maxV(2,:))

%% Calculate the Right COP V
COP = zeros(fraFp,2);
for i = 1:fraFp
    COP(i,1) = (-M2(i,2)+F2(i,1)*0.0178)/F2(i,3)+0.232; % COP coordination in global
    COP(i,2) = (M2(i,1)+F2(i,2)*0.0178)/F2(i,3)-0.254;
end

%% Transform the FP Frams to Vicon Frame V
footRd = F2(1:8:end,:);
COPfsVi = COP(1:8:end,:);
MfsVi = M2(1:8:end,:);
%% Calculate the Segment to Global Coordination Matrix V
[femSG,femGS,femDeg2,femDeg] = SegToGlo(Fem);
[tibSG,tibGS,tibDeg2,tibDeg] = SegToGlo(Tib);
[footSG,footGS,footDeg2,footDeg] = SegToGlo(Foot);

femDeg = FilSmoo(femDeg,2,fc,fsVi,'low');
tibDeg = FilSmoo(tibDeg,2,fc,fsVi,'low');
footDeg = FilSmoo(footDeg,2,fc,fsVi,'low');

femCOMg = zeros(fraVi,dim);
tibCOMg = zeros(fraVi,dim);
footCOMg = zeros(fraVi,dim);
for i = 1:fraVi
    femCOMg(i,:) = (0.567*Fem(i,:,4)+0.433*Fem(i,:,1));
    tibCOMg(i,:) = (0.567*Tib(i,:,4)+0.433*Tib(i,:,1));
    footCOMg(i,:) = (0.5*Foot(i,:,4)+0.5*Foot(i,:,1));
end

%% Calculate Segment Angular Velocity & Acceleration V
femV = vel_diff(femCOMg,1/fsVi); % GRS: linear velocity
tibV = vel_diff(tibCOMg,1/fsVi);
footV = vel_diff(footCOMg,1/fsVi);
femV = Fil(femV,2,fc,fsVi,'low');
tibV = Fil(tibV,2,fc,fsVi,'low');
footV = Fil(footV,2,fc,fsVi,'low');

femA = acc_diff(femCOMg,1/fsVi); % GRS: linear acceleration
tibA = acc_diff(tibCOMg,1/fsVi);
footA = acc_diff(footCOMg,1/fsVi);
femA = Fil(femA,2,fc,fsVi,'low');
tibA = Fil(tibA,2,fc,fsVi,'low');
footA = Fil(footA,2,fc,fsVi,'low');

femAngVs = SegAngVel(femDeg,1/fsVi).'; % LRS: angular velocity
tibAngVs = SegAngVel(tibDeg,1/fsVi).';
footAngVs = SegAngVel(footDeg,1/fsVi).';
femAngVs = Fil(femAngVs,2,fc,fsVi,'low');
tibAngVs = Fil(tibAngVs,2,fc,fsVi,'low');
footAngVs = Fil(footAngVs,2,fc,fsVi,'low');

femAngA = acc_diff(femDeg,1/fsVi); % GRS: angular acceleration 
tibAngA = acc_diff(tibDeg,1/fsVi);
footAngA = acc_diff(footDeg,1/fsVi);
femAngA = Fil(femAngA,2,fc,fsVi,'low');
tibAngA = Fil(tibAngA,2,fc,fsVi,'low');
footAngA = Fil(footAngA,2,fc,fsVi,'low');
for i = 1:fraVi % LRS: angular acceleration 
     femAngAs(i,:) = [cos(femDeg(i,2))*cos(femDeg(i,3)) sin(femDeg(i,3)) 0; -cos(femDeg(i,2))*sin(femDeg(i,3)) cos(femDeg(i,3)) 0;sin(femDeg(i,2)) 0 1]...
            *femAngA(i,:).';
     tibAngAs(i,:) = [cos(tibDeg(i,2))*cos(tibDeg(i,3)) sin(tibDeg(i,3)) 0; -cos(tibDeg(i,2))*sin(tibDeg(i,3)) cos(tibDeg(i,3)) 0;sin(tibDeg(i,2)) 0 1]...
             *tibAngA(i,:).'; 
     footAngAs(i,:) = [cos(footDeg(i,2))*cos(footDeg(i,3)) sin(footDeg(i,3)) 0; -cos(footDeg(i,2))*sin(footDeg(i,3)) cos(footDeg(i,3)) 0;sin(footDeg(i,2)) 0 1]...
            *tibAngA(i,:).';
end

%% Calculate Proximal Force in GRS and LRS
footRp(:,1) = footMass*footA(:,1)-footRd(:,1); 
footRp(:,2) = footMass*footA(:,2)-footRd(:,2);
footRp(:,3) = footMass*footA(:,3)+footMass*9.8-footRd(:,3);
footrd = GSRot(footGS,footRd);
footrp = GSRot(footGS,footRp);

tibRd = -footRp; % GRS: tib distal F = inverse foot proximal F
tibrd = GSRot(tibGS,tibRd); % LRS: tib distal F = inverse foot proximal F
tibRp(:,1) = tibMass*tibA(:,1)-tibRd(:,1); % Rp = ma - Rd
tibRp(:,2) = tibMass*tibA(:,2)-tibRd(:,2);
tibRp(:,3) = tibMass*tibA(:,3)+tibMass*9.8-tibRd(:,3);
tibrp = GSRot(tibGS,tibRp);

femRd = -tibRp;
femrd = GSRot(femGS,femRd);
femRp(:,1) = femMass*femA(:,1)-femRd(:,1); 
femRp(:,2) = femMass*femA(:,2)-femRd(:,2);
femRp(:,3) = femMass*femA(:,3)+femMass*9.8-femRd(:,3);
femrp = GSRot(femGS,femRp);

%% Calculate the Moment in GRS and LRS V
footmp = zeros(fraVi,dim);
tibmp = zeros(fraVi,dim);
femmp = zeros(fraVi,dim);

footmd = GSRot(footGS,MfsVi);
for i = 1:fraVi
    footmp(i,1) = footI(1)*footAngAs(i,1)-(footI(2)-footI(3))*footAngVs(i,2)*footAngVs(i,3)-footrd(i,2)*footL*0.5+footrp(i,2)*footL*0.5-footmd(i,1);
    footmp(i,2) = footI(2)*footAngAs(i,2)-(footI(3)-footI(1))*footAngVs(i,1)*footAngVs(i,3)+footrd(i,1)*footL*0.5-footrp(i,1)*footL*0.5-footmd(i,2);
    footmp(i,3) = footI(3)*footAngAs(i,3)-(footI(1)-footI(2))*footAngVs(i,1)*footAngVs(i,2)-footmd(i,3);
end
footMp = GSRot(footSG,footmp);
tibmd = GSRot(tibGS,-footMp);
for i = 1:fraVi
    tibmp(i,1) = tibI(1)*tibAngAs(i,1)-(tibI(2)-tibI(3))*tibAngVs(i,2)*tibAngVs(i,3)-tibrd(i,2)*tibL*0.567+tibrp(i,2)*tibL*0.433-tibmd(i,1);
    tibmp(i,2) = tibI(2)*tibAngAs(i,2)-(tibI(3)-tibI(1))*tibAngVs(i,1)*tibAngVs(i,3)+tibrd(i,1)*tibL*0.567-tibrp(i,1)*tibL*0.433-tibmd(i,2);
    tibmp(i,3) = tibI(3)*tibAngAs(i,3)-(tibI(1)-tibI(2))*tibAngVs(i,1)*tibAngVs(i,2)-tibmd(i,3);
end
tibMp = GSRot(tibSG,tibmp);
femmd = GSRot(femGS,-tibMp);
for i = 1:fraVi
    femmp(i,1) = femI(1)*femAngAs(i,1)-(femI(2)-femI(3))*femAngVs(i,2)*femAngVs(i,3)-femrd(i,2)*femL*0.567+femrp(i,2)*femL*0.433-femmd(i,1);
    femmp(i,2) = femI(2)*femAngAs(i,2)-(femI(3)-femI(1))*femAngVs(i,1)*femAngVs(i,3)+femrd(i,1)*femL*0.567-femrp(i,1)*femL*0.433-femmd(i,2);
    femmp(i,3) = femI(3)*femAngAs(i,3)-(femI(1)-femI(2))*femAngVs(i,1)*femAngVs(i,2)-femmd(i,3);
end
femMp = GSRot(femSG,femmp);

footmp = FilSmoo(footmp,2,fc,fsVi,'low')/BM;
tibmp = FilSmoo(tibmp,2,fc,fsVi,'low')/BM;
femmp = FilSmoo(femmp,2,fc,fsVi,'low')/BM;

footMp = FilSmoo(footMp,2,fc,fsVi,'low')/BM;
tibMp = FilSmoo(tibMp,2,fc,fsVi,'low')/BM;
femMp = FilSmoo(femMp,2,fc,fsVi,'low')/BM;

%% Calculate the Power
Vrad = deg2rad(V);
P = zeros(fraVi,dim);
P(:,1) = -femMp(:,2).*Vrad(:,1);
P(:,2) = tibMp(:,2).*Vrad(:,2);
P(:,3) = -footMp(:,2).*Vrad(:,3);

for i = 1:fraVi
    if i ~= 1 && i ~= fraVi
        W(i,:) = (P(i,:)+P(i-1,:))/2;
    else
        W(i,:) = (P(i,:)+P(i,:))/2;
    end
    W(i,:) = W(i,:)*1/fsVi;
end
W = FilSmoo(W,2,fc,fsVi,'low');
EccW = sum(W(eccStaVi:conStaVi-1,:),1);
ConW = sum(W(conStaVi:conEndVi,:),1);

WP = zeros(1,dim);
WN = zeros(1,dim);
for i = 1:dim
    for j = eccStaVi:conEndVi
        if W(j,i) >= 0
            WP(1,i) = WP(1,i) + W(j,i);
        else
            WN(1,i) = WN(1,i) + W(j,i);
        end
    end
end


%% Get the Peak Force, Impulse, Moment
maxFz = [max(tFz(eccSta:conSta-1)) max(tFz(conSta:conEnd))];
maxFz = (maxFz-BW)/BM;
netI = ([trapz(tFz(unSta:conSta-1)-BW)/fsFp trapz(tFz(conSta:conEnd)-BW)/fsFp])/BM;
Vto = (netI(2)-netI(1))/BM;
JHfp = 0.5*(Vto^2)/9.8;
maxM = [max(femMp(eccStaVi:conStaVi-1,2)) -min(tibMp(eccStaVi:conStaVi-1,2)) max(footMp(eccStaVi:conStaVi-1,2));max(femMp(conStaVi:conEndVi,2)) -min(tibMp(conStaVi:conEndVi,2)) max(footMp(conStaVi:conEndVi,2))];
maxM2 = [max(femmp(eccStaVi:conStaVi-1,2)) -min(tibmp(eccStaVi:conStaVi-1,2)) max(footmp(eccStaVi:conStaVi-1,2));max(femmp(conStaVi:conEndVi,2)) -min(tibmp(conStaVi:conEndVi,2)) max(footmp(conStaVi:conEndVi,2))];

if ROM(1) > ROM(2)
    Strategy(1) = 1;
else
    Strategy(1) = 2;
end
if maxM(2,1) > maxM(2,2)
    Strategy(2) = 1;
else
    Strategy(2) = 2;
end
if maxM2(2,1) > maxM2(2,2)
    Strategy(3) = 1;
else
    Strategy(3) = 2;
end
%% Print Kinetics Data
fprintf('Eccentric peak force=%.4f(N/kg)\t Concentric peak force=%.4f(N/kg)\n', maxFz)
fprintf('Eccentric net impulse=%.4f(N/kg)\t Concentric net impulse=%.4f(N/kg)\n', netI)
fprintf('Max. moment(eccentric):Hip=%.4f(Nm/kg)\t Knee=%.4f(Nm/kg)\t Ankle=%.4f(Nm/kg)\n', maxM(1,:))
fprintf('Max. moment(concentric):Hip=%.4f(Nm/kg)\t Knee=%.4f(Nm/kg)\t Ankle=%.4f(Nm/kg)\n', maxM(2,:))
fprintf('Max. moment(eccentric):Hip=%.4f(Nm/kg)\t Knee=%.4f(Nm/kg)\t Ankle=%.4f(Nm/kg)\n', maxM2(1,:))
fprintf('Max. moment(concentric):Hip=%.4f(Nm/kg)\t Knee=%.4f(Nm/kg)\t Ankle=%.4f(Nm/kg)\n', maxM2(2,:))
fprintf('Net work(eccentric):Hip=%.4f(Nm/kg)\t Knee=%.4f(Nm/kg)\t Ankle=%.4f(Nm/kg)\n', WN)
fprintf('Net work(concentric):Hip=%.4f(Nm/kg)\t Knee=%.4f(Nm/kg)\t Ankle=%.4f(Nm/kg)\n', WP)

%% Plot 
subplot(3,3,1); % Plot force, velCom, disCom and True COM with Each Phase
plot(tFz)
hold on
plot([eccSta eccSta],[-500 2500],"magenta")
plot([unSta unSta],[-500 2500])
plot([conSta conSta],[-500 2500],"magenta")
plot([conEnd conEnd],[-500 2500],"magenta")
title('Fz','fontweight','bold','fontsize',12)

subplot(3,3,2); 
plot(comVFp)
hold on
plot([eccSta eccSta],[-4 4],"magenta")
plot([unSta unSta],[-4 4])
plot([conSta conSta],[-4 4],"magenta")
plot([conEnd conEnd],[-4 4],"magenta")
title('Com Velocity','fontweight','bold','fontsize',12)

subplot(3,3,3); 
plot(Com)
hold on
plot([eccStaVi eccStaVi],[0 200],"magenta")
plot([conStaVi conStaVi],[0 200],"magenta")
plot([conEndVi conEndVi],[0 200],"magenta")
plot([conEndVi+maxComFra conEndVi+maxComFra],[0 200])
title('Com Displacement','fontweight','bold','fontsize',12)

subplot(3,3,4); 
plot(Ang(:,1),'red')
hold on
plot(Ang(:,2),'blue')
plot(Ang(:,3),'green')
plot([eccStaVi eccStaVi],[-100 100],"magenta")
plot([conStaVi conStaVi],[-100 100],"magenta")
plot([conEndVi conEndVi],[-100 100],"magenta")
title('Joint Angle','fontweight','bold','fontsize',12)

subplot(3,3,5); 
plot(V(:,1),'red')
hold on
plot(V(:,2),'blue')
plot(V(:,3),'green')
plot([eccStaVi eccStaVi],[-1000 1000],"magenta")
plot([conStaVi conStaVi],[-1000 1000],"magenta")
plot([conEndVi conEndVi],[-1000 1000],"magenta")
title('Joint Velocity','fontweight','bold','fontsize',12)

subplot(3,3,6); 
plot(A(:,1),'red')
hold on
plot(A(:,2),'blue')
plot(A(:,3),'green')
plot([eccStaVi eccStaVi],[-10000 20000],"magenta")
plot([conStaVi conStaVi],[-10000 20000],"magenta")
plot([conEndVi conEndVi],[-10000 20000],"magenta")
title('Joint Acceleration','fontweight','bold','fontsize',12)

subplot(3,3,7);
plot(femMp(:,2),'red')
hold on
plot(tibMp(:,2),'blue')
plot(footMp(:,2),'green')
plot([eccStaVi eccStaVi],[-5 5],"magenta")
plot([conStaVi conStaVi],[-5 5],"magenta")
plot([conEndVi conEndVi],[-5 5],"magenta")
title('Joint Moment','fontweight','bold','fontsize',12)

subplot(3,3,8);
plot(P(:,1),'red')
hold on
plot(P(:,2),'blue')
plot(P(:,3),'green')
plot([eccStaVi eccStaVi],[-20 20],"magenta")
plot([conStaVi conStaVi],[-20 20],"magenta")
plot([conEndVi conEndVi],[-20 20],"magenta")
title('Joint Power','fontweight','bold','fontsize',12)

subplot(3,3,9);
plot(W(:,1),'red')
hold on
plot(W(:,2),'blue')
plot(W(:,3),'green')
plot([eccStaVi eccStaVi],[-1 1],"magenta")
plot([conStaVi conStaVi],[-1 1],"magenta")
plot([conEndVi conEndVi],[-1 1],"magenta")
title('Joint Work','fontweight','bold','fontsize',12)


%cd(outputPath);

%T = table({ID},{Task},Strategy,JH,Depth,eccT,conT,maxFz,netI,maxAng,ROM,maxV(1,1:2),maxV(2,1:2),maxM(1,:),maxM(2,:),maxM2(1,:),maxM2(2,:),WN,WP);
%writetable(T,'Output Data.xlsx','WriteMode','append');

