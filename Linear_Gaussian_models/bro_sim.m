function simulated=bro_sim(t,gap)
% General Settings%
%t=0:0.01:100;%
step=t/gap;
int64(step);
A=eye(2);
B=eye(2);
C=cell(step,1);
for i=1:step
    oberr.c1(i)=1+0.001*randn(1);
    oberr.c2(i)=1+0.001*randn(1);
end
for k=1:step
    C{k,1}=[oberr.c1(k) 0 ;0 oberr.c2(k)];
end
d = 1.0e-4; % diameter in meters
eta = 1.002e-3; % viscosity of water in SI units (Pascal-seconds)
kB = 1.38e-23; % Boltzmann constant
T = 300; % Temperature in degrees Kelvin
D = kB * T / (3 * pi * eta * d);
k=sqrt(2*D*gap);
par=brown(k,step);
x1=zeros(step,1);
x2=zeros(step,1);
y.m_backward=cell(1,step);
x.m_forward=cell(1,step);
x.V_forward=cell(1,step);
x.W_forward=cell(1,step);
x.m_p_forward=cell(1,step);
x.V_p_forward=cell(1,step);
Z.m_forward=cell(1,step);
Z.V_forward=cell(1,step);
U.m_forward=[0;0];
U.V_forward=eye(2);
x.m_backward=cell(1,step);
x.W_backward=cell(1,step);
x.Wm_backward=cell(1,step);
x.m_p_backward=cell(1,step);
x.W_p_backward=cell(1,step);
x.Wm_p_backward=cell(1,step);
x.m_dp_backward=cell(1,step);
x.W_dp_backward=cell(1,step);
x.Wm_dp_backward=cell(1,step);
Z.m_backward=cell(1,step);
Z.W_backward=cell(1,step);
U.W_forward=inv(U.V_forward);
y.W_backward=inv(10^-10*eye(2));
% Generating observation %
for i=1:step
    par.x1(i)=10^4*par.x1(i);
    par.x2(i)=10^4*par.x2(i);
end
for i=1:step
    y.m_backward{1,i}=[par.x1(i);par.x2(i)];
end
% Forward Recursion
x.V_forward{1,1}=10^-10*eye(2);
x.W_forward{1,1}=inv(x.V_forward{1,1});
x.m_forward{1,1}=y.m_backward{1,1};
for i=2:step
    G=inv(10^-10*eye(2)+C{i-1,1}*x.V_forward{1,i-1}*C{i-1,1}');
    x.V_p_forward{1,i-1}=x.V_forward{1,i-1}-x.V_forward{1,i-1}*C{i-1,1}'*G*C{i-1}*x.V_forward{1,i-1};
    x.m_p_forward{1,i-1}=x.m_forward{1,i-1}+x.V_forward{1,i-1}*C{i-1,1}'*G*(y.m_backward{1,i-1}-C{i-1,1}*x.m_forward{1,i-1});
    Z.m_forward{1,i}=A*x.m_p_forward{1,i-1};
    Z.V_forward{1,i}=A*x.V_p_forward{1,i-1}*A';
    x.m_forward{1,i}=Z.m_forward{1,i}+B*U.m_forward;
    x.V_forward{1,i}=Z.V_forward{1,i}+B*U.V_forward*B';
    x.W_forward{1,i}=inv(x.V_forward{1,i});
end
% Backward Recursion %
x.m_p_backward{1,step}=[0;0];
x.W_p_backward{1,step}=inv(10^8*eye(2));
x.Wm_p_backward{1,step}=x.W_p_backward{1,step}*x.m_p_backward{1,step};
for i=step:-1:1
    x.Wm_dp_backward{1,i}=C{i,1}'*y.W_backward*y.m_backward{1,i};
    x.W_dp_backward{1,i}=C{i,1}'*y.W_backward*C{i,1};
    x.Wm_backward{1,i}=x.Wm_p_backward{1,i}+x.Wm_dp_backward{1,i};
    x.W_backward{1,i}=x.W_p_backward{1,i}+x.W_dp_backward{1,i};
    x.m_backward{1,i}=x.W_backward{1,i}\x.Wm_backward{1,i};
    Z.m_backward{1,i}=x.m_backward{1,i}+B*U.m_forward;
    H=inv(U.W_forward+B'*x.W_backward{1,i}*B);
    Z.W_backward{1,i}=x.W_backward{1,i}-x.W_backward{1,i}*B*H*B'*x.W_backward{1,i};
    if i==1
        x.Wm_backward{1,i}=x.Wm_p_backward{1,i}+x.Wm_dp_backward{1,i};
        break;
    end
    x.W_p_backward{1,i-1}=A'*Z.W_backward{1,i}*A;
    x.Wm_p_backward{1,i-1}=A'*Z.W_backward{1,i}*Z.m_backward{1,i};
end
% Error terms %
err=zeros(step,1);
err1=zeros(step,1);
err2=zeros(step,1);
parfor i=1:step
    result=particle(i,x);
    x1(i,1)=result.x1;
    x2(i,1)=result.x2;
    err(i,1)=sqrt((par.x1(i,1)-x1(i,1)).^2+(par.x2(i,1)-x2(i,1)).^2);
    err1(i,1)=(par.x1(i,1)-x1(i,1));
    err2(i,1)=par.x2(i,1)-x2(i,1);
end
simulated=struct('x1',x1,'x2',x2);
simulated.RMSE=sqrt(mean((par.x1-simulated.x1).^2+(par.x2-simulated.x2).^2));
subplot(3,2,[1,2]);
histogram(err);
ylabel('Number of Error Terms');
xlabel('Euclidean Distance between Generated and Estimated position (mm)');
title([num2str(t), 's with ', num2str(gap), 's time gap Brownian motion simulation']);
subplot(3,2,3);
histogram(err1);
ylabel('Number of Error Terms');
xlabel('x - estimated x (mm)');
subplot(3,2,4);
histogram(err2);
ylabel('Number of Error Terms');
xlabel('y - estimated y (mm)');
subplot(3,2,5);
plot(par.x1,par.x2)
ylabel('Y Position');
xlabel('X Position');
title('Generated position versus time in 2D (mm)');
subplot(3,2,6);
plot(simulated.x1,simulated.x2);
ylabel('Y Position');
xlabel('X Position');
title('Estimated position versus time in 2D (mm)');
end