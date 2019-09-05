function tra = motionsim(mass,init_x1,init_x2,init_vx1,init_vx2,time,delta,Area)
% Generate a specified motion %
tra=struct();
step=time/delta;
int64(step);
tra.x1=zeros(step,1);
tra.x2=zeros(step,1);
tra.vx1=zeros(step,1);
tra.vx2=zeros(step,1);
tra.ax1=zeros(step,1);
tra.ax2=zeros(step,1);
tra.x1(1)=init_x1;
tra.x2(1)=init_x2;
g=9.81;%m/s^2 (acceleration due to gravity)
C=.5; %Drag Coefficient of a sphere
rho= 1.2; %kg/m^3 (density of air)
tra.vx1(1)=init_vx1;
tra.vx2(1)=init_vx2;
for i=1:step-1
    tra.ax1(i)=0-(C*rho/2*(tra.vx1(i))^2*Area)/mass;
    if tra.vx2(i)>0
        tra.ax2(i)=-g-(C*rho/2*(tra.vx2(i))^2*Area)/mass;
    else
        tra.ax2(i)=-g+(C*rho/2*(tra.vx2(i))^2*Area)/mass;
    end
    if(tra.x2(i)<-0.0001)
        s=(tra.x2(i)-tra.x2(i-1))/(tra.x1(i)-tra.x1(i-1));
        m=tra.x2(i)-s*tra.x1(i);
        tra.x1(i)=-m/s;
        tra.x2(i)=0;
        for k=i:step-1
            tra.x1(k+1)=tra.x1(k);
            tra.vx1(k+1)=0;
            tra.vx2(k+1)=0;
            tra.ax1(k+1)=0;
            tra.ax1(k+1)=0;
        end
        break;
    end
    tra.vx1(i+1)=tra.vx1(i)+tra.ax1(i)*delta;
    tra.vx2(i+1)=tra.vx2(i)+tra.ax2(i)*delta;
    tra.x1(i+1)=tra.x1(i)+tra.vx1(i)*delta;%+0.5*a.x1*delta^2;
    tra.x2(i+1)=tra.x2(i)+tra.vx2(i)*delta;%+0.5*a.x2*delta^2;
end
plot(tra.x1,tra.x2);
end