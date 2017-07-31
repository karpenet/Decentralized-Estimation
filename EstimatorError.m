%KarpeDiem
%Developed by Kedar Prasad Karpe

%This code simulates the tracking error of the dependent nodes thet hgd
%in a decentralized network system based on the decentralized state 
%estimator developed by Francesca Boem, Lorenzo Sabattini and Cristian 
%Secchi 

%Dependencies: topologies.m

clear all

%% Define Meta Variables

duration = 25;  %Duration of simulation
dt = 0.1;   %Step size

global trajectory_scale;
trajectory_scale = .001;

[I,W,Lv,n,T,J,J0] = topologies(3);    %Dependencies from set topologies

MaxEigenValue = 0;  %For pole placement of estimator dynamics

%% Define Network Matrices

NumberAgentScalar=size(I,1);
WeightsMatrix=diag(W);
LaplacianMatrix=I*WeightsMatrix*transpose(I);

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
C=-Tfl'*LaplacianMatrix*Tfl;

Nd = size(A,1); %Number of Dependent Robots
Ni = size(B,2); %Number of Independent Robots

%% Generate Estimator Gain Dependencies

%Generate matrix F
[F] = -lqr(A,B,eye(size(A,1)),eye(size(B,2)))';


%Generate Fe Matrix
Fe = kron(eye(Ni),F');


%Generate FeTilda Matrix
FeTilda = [];
FeTildaRow = [];

for index_J1 = 1:Ni
    Ftemp = F';
    Ftemp(index_J1,:)=0;
    TildaBlock = F'-Ftemp;
    FeTildaRow = [FeTildaRow TildaBlock];
end

for index_I1 = 1:Ni
    FeTilda = [FeTilda; FeTildaRow];
end


%% Generate estimator gain matrix Ke

Agen = A + B*F';

GainCond = -norm(kron(eye(Ni),B) * FeTilda);    %Estimator Gain Condition LambdaIMin

eigen = GainCond + ((rand(1,size(A,1))).*(MaxEigenValue-GainCond));   %Generate Eigen Values from Condition


Ki = [];
for indexKi = 1:Ni
    ki = place(Agen', B(:,indexKi), eigen)';
    Ki = [Ki ki];
end
 

%% Error Estimator Dynamics Dependencies

 Ae = kron(eye(Ni), A);
 Be = kron(eye(Ni), B);
 
 ExoMatrixTemp = [];
 for LambdaIndex=1:Ni
    MiniJordanBlock= A + Ki(:, LambdaIndex) * B(:,LambdaIndex)';
    ExoMatrixTemp = blkdiag(ExoMatrixTemp,MiniJordanBlock);
end;

Lambda = ExoMatrixTemp;

IeTilda = [];
IeTildaRow = [];

for index_J2 = 1:Ni
    TildaBlock = ones(Ni);
    TildaBlock(index_J2, index_J2) = 0;
    IeTildaRow = [IeTildaRow TildaBlock];
end

for index_I2 = 1:Ni
    IeTilda = [IeTilda; IeTildaRow];
end

BTilda = Be*(eye(size(IeTilda,1),size(IeTilda,2))-IeTilda);

Psi = Lambda + kron(eye(Ni),B*F');

%% Estimator Error Dynamics

AErrTemp = Lambda + BTilda*Fe;
AErr = (Psi - (kron(eye(Ni),B)*FeTilda));
sys = ss(AErr, zeros(size(AErr,1),size(AErr,2)),eye(size(AErr,1)),0);
[~, t, EstimatorError] = initial(sys,round(rand(1,size(AErr,1))), 0:dt:duration);

%% Plot Estimation Error
plot(t,EstimatorError)
hold on
axis square
