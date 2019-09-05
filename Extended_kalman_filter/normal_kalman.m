function normal_kalman(gap,time)
    rng('shuffle');
    step=time/gap;
    int64(step);
    d=zeros(1,step);
    for i=1:step
        d(i)=i;
    end
    var_a=0.5;
    var_z=0.2;
    F=[1 gap;0 1];
    G=[0.5*gap^2;gap];
    Q=[0.25*gap^4 0.5*gap^3;0.5*gap^3 gap^2]*var_a;
    H=[1 0];
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
    for i=2:step
        truth{1,i}=F*truth{1,i-1}+G*sqrt(var_a)*randn(1);%true track
        x.vt(i)=truth{1,i}(2,1);
        x.lt(i)=truth{1,i}(1,1);
    end
    for i=1:step
        z{1,i}=H*truth{1,i}+sqrt(var_z)*randn(1);%generated observation
        x.lo(i)=z{1,i}(1,1);
    end
    x.posteriori{1,1}=truth{1,1}; %initial state
    P.posteriori{1,1}=zeros(2); %know exactly initial state -> 0 covariance
    err_l=zeros(1,step);
    err_v=zeros(1,step);
    for i=2:step
        x.priori{1,i}=F*x.posteriori{1,i-1};
        P.priori{1,i}=F*P.posteriori{1,i-1}*F.'+Q;
        y{1,i}=z{1,i}-H*x.priori{1,i};
        S{1,i}=H*P.priori{1,i}*H.'+R;
        K{1,i}=P.priori{1,i}*H.'*inv(S{1,i});
        x.posteriori{1,i}=x.priori{1,i}+K{1,i}*y{1,i};
        P.posteriori{1,i}=(eye(2)-K{1,i}*H)*P.priori{1,i};
        x.lp(i)=x.posteriori{1,i}(1,1);
        x.vp(i)=x.posteriori{1,i}(2,1);
        err_l(i)=abs(x.lt(i)-x.lp(i));
        err_v(i)=abs(x.vt(i)-x.vp(i));
    end
    subplot(2,2,1);
    plot(x.lt,d);
    hold on;
    plot(x.lo,d);
    hold on;
    plot(x.lp,d);
    hold off;
    subplot(2,2,2);
    plot(x.vt,d);
    hold on;
    plot(x.vp,d);
    subplot(2,2,3);
    histogram(err_l);
    subplot(2,2,4);
    histogram(err_v);
    
end