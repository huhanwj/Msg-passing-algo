function tra = unknown_m(mass,init_x1,init_x2,init_vx1,init_vx2,time,delta,Area)
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
    C=.5; %Drag Coefficient of a sphere
    rho= 1.2; %kg/m^3 (density of air)
    tra.vx1(1)=init_vx1;
    tra.vx2(1)=init_vx2;
    for i=1:step-1
        if i<step/4
            if tra.vx1(i)>0
                tra.ax1(i)=(-(C*rho/2*(tra.vx1(i))^2*Area))/mass;
            else
                tra.ax1(i)=((C*rho/2*(tra.vx1(i))^2*Area))/mass;
            end
            if tra.vx2(i)>0
                tra.ax2(i)=(-(C*rho/2*(tra.vx2(i))^2*Area))/mass;
            else
                tra.ax2(i)=((C*rho/2*(tra.vx2(i))^2*Area))/mass;
            end
        elseif i>=step/4 && i<step/2
            if tra.vx1(i)>0
                tra.ax1(i)=(log(3)*mass-(C*rho/2*(tra.vx1(i))^2*Area))/mass;
            else
                tra.ax1(i)=(log(3)*mass+(C*rho/2*(tra.vx1(i))^2*Area))/mass;
            end
            if tra.vx2(i)>0
                tra.ax2(i)=(log(2)*mass-(C*rho/2*(tra.vx2(i))^2*Area))/mass;
            else
                tra.ax2(i)=(log(2)*mass+(C*rho/2*(tra.vx2(i))^2*Area))/mass;
            end
        elseif i>=step/2 && i<4*step/5
            if tra.vx1(i)>0
                tra.ax1(i)=(-0.5*mass-(C*rho/2*(tra.vx1(i))^2*Area))/mass;
            else
                tra.ax1(i)=(-log(5)*mass+(C*rho/2*(tra.vx1(i))^2*Area))/mass;
            end
            if tra.vx2(i)>0
                tra.ax2(i)=(-0.3*mass-(C*rho/2*(tra.vx2(i))^2*Area))/mass;
            else
                tra.ax2(i)=(-log(1.5)*mass+(C*rho/2*(tra.vx2(i))^2*Area))/mass;
            end
        else
            if tra.vx1(i)>0
                tra.ax1(i)=(mass-(C*rho/2*(tra.vx1(i))^2*Area))/mass;
            else
                tra.ax1(i)=(mass+(C*rho/2*(tra.vx1(i))^2*Area))/mass;
            end
            if tra.vx2(i)>0
                tra.ax2(i)=(mass-(C*rho/2*(tra.vx2(i))^2*Area))/mass;
            else
                tra.ax2(i)=(mass+(C*rho/2*(tra.vx2(i))^2*Area))/mass;
            end
        end
        %if(tra.x2(i)<-0.0001)
        %    s=(tra.x2(i)-tra.x2(i-1))/(tra.x1(i)-tra.x1(i-1));
        %    m=tra.x2(i)-s*tra.x1(i);
        %    tra.x1(i)=-m/s;
        %    for k=i:step-1
        %        tra.x1(k+1)=tra.x1(k);    
        %    end
        %    break;
        %end
        tra.vx1(i+1)=tra.vx1(i)+tra.ax1(i)*delta;
        tra.vx2(i+1)=tra.vx2(i)+tra.ax2(i)*delta;
        tra.x1(i+1)=tra.x1(i)+tra.vx1(i)*delta;%+0.5*a.x1*delta^2;
        tra.x2(i+1)=tra.x2(i)+tra.vx2(i)*delta;%+0.5*a.x2*delta^2;
    end
    %plot(tra.x1,tra.x2);
    end