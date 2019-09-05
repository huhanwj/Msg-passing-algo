function simulated=motion_est(mass,init_x1,init_x2,init_vx1,init_vx2,time,delta,Area)
step=time/delta;
int64(step);
rng('shuffle');
% variable setup %
x1=zeros(1,step);
vx1=zeros(1,step);
x2=zeros(1,step);
vx2=zeros(1,step);
% A B C Matrix setup %
A=[1 delta 0 0;0 1 0 0;0 0 1 delta;0 0 0 1];
B=[0 0 0 0;0 delta 0 0;0 0 0 0;0 0 0 delta];
C=cell(step,1);
for i=1:step
    oberr.c1(i)=1;
    oberr.c2(i)=1;
    oberr.c3(i)=1;
    oberr.c4(i)=1;
end
for k=1:step
    C{k}=[oberr.c1(k) 0 0 0;0 oberr.c2(k) 0 0;0 0 oberr.c3(k) 0;0 0 0 oberr.c4(k)];
end
% Generating comparison data %
par=motionsim(mass,init_x1,init_x2,init_vx1,init_vx2,time,delta,Area);
%par=unknown_m(mass,init_x1,init_x2,init_vx1,init_vx2,time,delta,Area);
% Generating observation %
y.m_backward=cell(step,1);
y.V_backward=10^(-3)*eye(4);
for i=1:step
    y.m_backward{i,1}=[par.x1(i);par.vx1(i);par.x2(i);par.vx2(i)];
end
% Generate input acceleration %
U.m_forward=cell(1,step);
U.V_forward=10^(-20)*eye(4);
for i=1:step
    U.m_forward{1,i}=[0;par.ax1(i);0;par.ax2(i)];
end
% Forward Recursion %
% message variables Setup %
x.m_forward=cell(1,step);
x.m_forward=cell(1,step);
x.V_forward=cell(1,step);
x.W_forward=cell(1,step);
x.m_p_forward=cell(1,step);
x.V_p_forward=cell(1,step);
Z.m_forward=cell(1,step);
Z.V_forward=cell(1,step);
% Initial condition %
x.V_forward{1,1}=10^-20*eye(4);
x.W_forward{1,1}=inv(x.V_forward{1,1});
x.m_forward{1,1}=[init_x1 ; init_vx1 ; init_x2 ; init_vx2];
% message calculation %
for i=2:step
    G=inv(y.V_backward+C{i-1,1}*x.V_forward{1,i-1}*C{i-1,1}');
    x.V_p_forward{1,i-1}=x.V_forward{1,i-1}-x.V_forward{1,i-1}*C{i-1,1}'*G*C{i-1,1}*x.V_forward{1,i-1};
    x.m_p_forward{1,i-1}=x.m_forward{1,i-1}+x.V_forward{1,i-1}*C{i-1,1}'*G*(y.m_backward{i-1,1}-C{i-1,1}*x.m_forward{1,i-1});
    Z.m_forward{1,i}=A*x.m_p_forward{1,i-1};
    Z.V_forward{1,i}=A*x.V_p_forward{1,i-1}*A';
    x.m_forward{1,i}=Z.m_forward{1,i}+B*U.m_forward{1,i};
    x.V_forward{1,i}=Z.V_forward{1,i}+B*U.V_forward*B';
    x.W_forward{1,i}=inv(x.V_forward{1,i});
end
% Backward Recursion %
% message variables setup %
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
y.W_backward=inv(y.V_backward);
% Initial Condition %
x.m_p_backward{1,step}=[0;0;0;0];
x.W_p_backward{1,step}=inv(10^20*eye(4));
x.Wm_p_backward{1,step}=x.W_p_backward{1,step}*x.m_p_backward{1,step};
% message calculation %
for i=step:-1:1
    x.Wm_dp_backward{1,i}=C{i,1}'*y.W_backward*y.m_backward{i,1};
    x.W_dp_backward{1,i}=C{i,1}'*y.W_backward*C{i,1};
    x.Wm_backward{1,i}=x.Wm_p_backward{1,i}+x.Wm_dp_backward{1,i};
    x.W_backward{1,i}=x.W_p_backward{1,i}+x.W_dp_backward{1,i};
    x.m_backward{1,i}=x.W_backward{1,i}\x.Wm_backward{1,i};
    Z.m_backward{1,i}=x.m_backward{1,i}+B*U.m_forward{1,i};
    H=inv(U.W_forward+B'*x.W_backward{1,i}*B);
    Z.W_backward{1,i}=x.W_backward{1,i}-x.W_backward{1,i}*B*H*B'*x.W_backward{1,i};
    if i==1
        break;
    end
    x.W_p_backward{1,i-1}=A'*Z.W_backward{1,i}*A;
    x.Wm_p_backward{1,i-1}=A'*Z.W_backward{1,i}*Z.m_backward{1,i};
end
% Error Terms setup %
err=zeros(step,1);
err1=zeros(step,1);
err2=zeros(step,1);
parfor i=1:step
    output=parsim(i,x);
    x1(1,i)=output.x1;
    vx1(1,i)=output.vx1;
    x2(1,i)=output.x2;
    vx2(1,i)=output.vx2;
    err(i,1)=sqrt((par.x1(i,1)-x1(1,i)).^2+(par.x2(i,1)-x2(1,i)).^2)
    err1(i,1)=par.x1(i,1)-x1(1,i);
    err2(i,1)=par.x2(i,1)-x2(1,i);
end
simulated=struct('x1',x1,'x2',x2,'err',err,'vx1',vx1,'vx2',vx2);
% simulated.RMSE=sqrt(mean((par.x1-simulated.x1)^2+(par.x2-simulated.x2)^2));
subplot(3,2,[1,2]);
%edges=zeros(1,100);
%edges(1)=0;
%for i=2:100
%    edges(i)=edges(i-1)+0.06;
%end
%edges = [0:0.0002:0.002];
histogram(simulated.err);
ylabel('Error Terms');
xlabel('Euclidean Distance between Generated and Estimated position (m)');
subplot(3,2,3);
histogram(err1);
ylabel('Number of Error Terms');
xlabel('x - estimated x (m)');
subplot(3,2,4);
histogram(err2);
ylabel('Number of Error Terms');
xlabel('y - estimated y (m)');
subplot(3,2,5);
plot(par.x1,par.x2)
ylabel('Y Position');
xlabel('X Position');
title('Generated position versus time in 2D (m)');
subplot(3,2,6);
plot(simulated.x1,simulated.x2);
ylabel('Y Position');
xlabel('X Position');
title('Estimated position versus time in 2D (m)');
end
