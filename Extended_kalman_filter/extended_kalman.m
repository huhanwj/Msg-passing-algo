function extended_kalman(gap,time)
    rng('shuffle');
    step=time/gap;
    int64(step);
    d=zeros(1,step);
    for i=1:step
        d(i)=i;
    end
    var_a=0.1;
    var_z=0.001;
    a=1;
    b=5*10^(-9);
    f=cell(1,step);
    F=cell(1,step);
    h=cell(1,step);
    H=cell(1,step);
    G=[0.5*gap^2;gap];
    Q=[0.25*gap^4 0.5*gap^3;0.5*gap^3 gap^2]*var_a;
    R=[var_z];
    x.priori=cell(1,step);
    x.posteriori=cell(1,step);
    x.vt=zeros(1,step);
    x.lt=zeros(1,step);
    x.vp=zeros(1,step);
    x.lp=zeros(1,step);
    x.lo=zeros(1,step);
    P.priori=cell(1,step);
    P.posteriori=cell(1,step);
    y=cell(1,step);
    S=cell(1,step);
    K=cell(1,step);
    truth=cell(1,step);
    truth{1,1}=[0;0];
    x.vt(1)=0;
    x.lt(1)=0;
    x.lo(1)=0;
    x.vp(1)=0;
    x.lp(1)=0;
    z=cell(1,step);
    f{1,1}=[1 gap;0 1];
    for i=2:step
        f{1,i}=[1 gap;0 1+gap*a*truth{1,i-1}(2,1)^2];
        truth{1,i}=f{1,i}*truth{1,i-1}+G*sqrt(var_a)*randn(1);%true track
        x.vt(i)=truth{1,i}(2,1);
        x.lt(i)=truth{1,i}(1,1);
    end
    for i=1:step
        h{1,i}=[1 b*(truth{1,i}(2,1))];
        z{1,i}=h{1,i}*truth{1,i}+sqrt(var_z)*randn(1);%generated observation
        x.lo(i)=z{1,i}(1,1);
    end
    x.posteriori{1,1}=truth{1,1}; %initial state
    P.posteriori{1,1}=zeros(2); %know exactly initial state -> 0 covariance
    err_l=zeros(1,step);
    err_v=zeros(1,step);
    for i=2:step
        x.priori{1,i}=f{1,i}*x.posteriori{1,i-1};
        F{1,i}=[1 gap;0 1+0.03*gap*x.posteriori{1,i-1}(2,1)^2];
        P.priori{1,i}=F{1,i}*P.posteriori{1,i-1}*F{1,i}.'+Q;
        y{1,i}=z{1,i}-h{1,i}*x.priori{1,i};
        H{1,i}=[1 2*b*x.priori{1,i}(2,1)];
        S{1,i}=H{1,i}*P.priori{1,i}*H{1,i}.'+R;
        K{1,i}=P.priori{1,i}*H{1,i}.'/S{1,i};
        x.posteriori{1,i}=x.priori{1,i}+K{1,i}*y{1,i};
        P.posteriori{1,i}=(eye(2)-K{1,i}*H{1,i})*P.priori{1,i};
        x.lp(i)=x.posteriori{1,i}(1,1);
        x.vp(i)=x.posteriori{1,i}(2,1);
        err_l(i)=abs(x.lt(i)-x.lp(i));
        err_v(i)=abs(x.vt(i)-x.vp(i));
    end
    subplot(2,2,1);
    plot(d,x.lt);
    hold on;
    plot(d,x.lo);
    hold on;
    plot(d,x.lp);
    hold off;
    xlabel({'Time T',[num2str(step),' ',num2str(gap),'s-step']});
    ylabel('Position (m)');
    legend("Actual Pos","Observation","Predicted Pos");
    legend('boxoff');
    subplot(2,2,2);
    plot(d,x.vt);
    hold on;
    plot(d,x.vp);
    hold off;
    xlabel({'Time T',[num2str(step),' ',num2str(gap),'s-step']});
    ylabel('Velocity (m\s)');
    legend("Actual Velocity","Predicted Velocity");
    legend('boxoff');
    subplot(2,2,3);
    histogram(err_l);
    xlabel('Difference between actual and predicted position');
    ylabel('Number of error terms');
    subplot(2,2,4);
    histogram(err_v);
    xlabel('Difference between actual and predicted velocity');
    ylabel('Number of error terms');
end