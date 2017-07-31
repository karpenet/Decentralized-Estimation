%KarpeDiem
%Developed by Kedar Prasad Karpe

%This code simulates the tracking error of the dependent nodes
%in a decentralized network system based on the decentralized state 
%estimator developed by Francesca Boem, Lorenzo Sabattini and Cristian 
%Secchi 

%Dependencies: slvRegEqn.m, compMotMat.m, compChiMat.m, compWMat.m,
%compHMat.m, topologies.m

clear all
%% Define Meta Variables
duration = 25;  %Duration of simulation
dt = 0.1;   %Step size

global trajectory_scale;
trajectory_scale = .001;

[I,W,Lv,n,T,J,J0] = topologies(3);
%% Define Network Matrices

NumberAgentScalar=size(I,1);
WeightsMatrix=diag(W);
LaplacianMatrix=I*WeightsMatrix*transpose(I);
plot(graph(LaplacianMatrix))

In=eye(NumberAgentScalar);
Pf=[];
Tfl=[];

for m=1:NumberAgentScalar
    if m~=Lv
        Pf=[Pf,In(:,m)];
    else
        Tfl=[Tfl,In(:,m)];
    end;
end;

%% Define Nodes

A=-Pf'*LaplacianMatrix*Pf;              
B=-Pf'*LaplacianMatrix*Tfl;

Nd = size(A,1); %Number of Dependent Robots
Ni = size(B,2); %Number of Independent Robots
%% Other variables
poles = -20*rand(1,Nd); %Poles matrix for pole placement
n_har_d = 2*n + 1;
GCtrl = [];

%% Setpoint Harmonics
harfunc_init = zeros(1,n);


%Generate Initial Condition Matrix for setpoint Diff Eqn
for x = 1:n_har_d
    if rem(x,2)==0
        harfunc_init(x)=0;
    else
        harfunc_init(x)=1;
    end
end
harfunc_init = transpose(harfunc_init);


%Generate Block Diagonal Matrix G

ExoVectorTemp=[1];
ExoMatrixTemp=0;

for FirstIndexScalar=1:fix(n)
    MiniJordanBlock=[0,FirstIndexScalar*2*pi/T;-FirstIndexScalar*2*pi/T,0];
    ExoMatrixTemp = blkdiag(ExoMatrixTemp,MiniJordanBlock);
    ExoVectorTemp=[ExoVectorTemp;[0;1]];
end;

GCtrl=ExoMatrixTemp;


%% Calculate Input Dependencies
[F] = -lqr(A,B,eye(size(A,1)),eye(size(B,2)))';

[Pi,Gamma]=slvRegEqs(A,B,GCtrl,J,J0);

%% Generate Harmonic Equation

Ag = [A + B * F', -B * F', B * (Gamma - F' * Pi);...
        zeros(size(A, 1), size(B * F', 2)), A + F * B', zeros(size(A, 1), size(B * (Gamma - F' * Pi), 2));... 
        zeros(size(GCtrl, 1), size(A + B * F', 2)), zeros(size(GCtrl, 1), size(B * F', 2)), GCtrl];

sys=ss(Ag,zeros(size(Ag,1),size(Ag,2)),eye(size(Ag,1)),zeros(size(Ag,1),size(Ag,2)));

[y,t,response]=initial(sys,[max(max(J)) * rand(2*size(A,1),1);harfunc_init],0:dt:duration);
        
%% dynamics

setdyn=response(:,1:size(A,1));
hate=response(:,size(A,1)+1:2*size(A,1));
harfunc=response(:,2*size(A,1)+1:size(response,2));
u=(F' * (setdyn' - hate') + (Gamma - F' * Pi) * harfunc')';
depdyn=harfunc*Pi';

errmat=depdyn-setdyn;

figure(2)
plot(t,errmat)
hold on
axis square
