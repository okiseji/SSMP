global WINDOW_SIZE MAXITER LENGTH NUMBER TOL;
WINDOW_SIZE=50;
MAXITER=1000;
TOL=10e-7;

data = load(fullfile('./DataSets',  'FF25EU.mat'));
data=data.data;
plot(data);
[LENGTH,NUMBER]=size(data);
rou=0.066;
WEALTH=1000;
We=[0];
p=1;
% OMV
for i=WINDOW_SIZE+1:LENGTH
    AP=data(i-WINDOW_SIZE:i-1,:);
    [mu,sigma]=ER(AP);
    A=[mu.';ones(NUMBER,1).'];
    b=[rou;1];
    D=[A;-A;eye(NUMBER)];
    d=[b;-b;zeros(NUMBER,1)];
    w=zeros(NUMBER,1);
    w(1,1)=1;
    y=zeros(NUMBER+4,1);
    y(2,1)=1;
    L=2*norm(sigma);
    beta_upperBound=2/L;
    beta=beta_upperBound/2;
    omega_upperBound=(2*(2-beta*L))/(4*beta*norm(D)*norm(D)+L*(2-beta*L));
    omega=omega_upperBound/2;
    for l=1:MAXITER
        nw=w-beta*(2*sigma*w+D.'*y);
        tc=(1/omega)*y+D*(2*nw-w);
        ny=omega*(tc-max(tc,d));
        w=nw;
        y=ny;
    end
    pv=data(i-1,:).';
    nv=data(i,:).';
    num=(WEALTH*w)./pv;
    WEALTH=nv.'*num;
    We(p)=WEALTH;
    p=p+1;
end

plot(We);
% average return and return covariance
function [mu,sigma]=ER(AP)
    global WINDOW_SIZE
    PV=AP(2:WINDOW_SIZE,:);
    FV=AP(1:WINDOW_SIZE-1,:);
    RM=(PV-FV)./FV;
    mu=mean(RM).';
    sigma=cov(RM);
end
