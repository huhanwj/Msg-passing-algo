rng('shuffle');   %random number seed
map=zeros(1,100);
x.sim=zeros(1,100);
result=zeros(1,100);
for i=1:100
    x_msg1=[1.5 2]; %Variable X
    z1_msg=[0 1];   %AWGN Z_1
    z2_msg=[0 1.5]; %AWGN Z_2
    x.sim(1,i)=x_msg1(1,1)+sqrt(x_msg1(1,2))*rand(1); %X vector generated from the setting
    y1=x_msg1(1,1)+sqrt(z1_msg(1,2))*rand(1);   %simulated observed Y_1
    y2=x_msg1(1,1)+sqrt(z2_msg(1,2))*rand(1);   %simluated observed Y_2
    y1_msg=[y1 0];  %message passed from y_1 to +
    y2_msg=[y2 0];  %message passed from y_2 to +
    x1_msg=z1_msg+y1_msg;   %message passed from + to = L.H.S.
    x2_msg=z2_msg+y2_msg;   %message passed from + to = R.H.S.
    W_x1=1/x_msg1(1,2); 
    %Weight matrix of in direction message of X (1x1 matrix since X is a scalar)
    W_x2=1/z1_msg(1,2)+1/z2_msg(1,2);   
    %Weight matrix of opposite direction message of X (1x1 matrix)
    W_x=W_x1+W_x2;  %Weight matrix of X
    V_x=1/W_x;      %Covariance matrix of X
    wm_1=W_x1*x_msg1(1,1);  %Wxmx matrix of in direction message of X
    wm_2=1/z1_msg(1,2)+x1_msg(1,1)+1/z2_msg(1,2)+x2_msg(1,1);   
    %Wxmx matrix of opposite direction message of X
    m_x=V_x*wm_1+V_x*wm_2; %mean of X
    x_marginal=@(x)-normpdf(x,m_x,V_x); %posterior probability of X given y1,y2
    map(1,i)=fminbnd(x_marginal,-10,10);%find MAP
    result(1,i)=sqrt((map(1,i)-x.sim(1,i))^2);
end
RMSE=sqrt(mean((map-x.sim).^2))   %RMSE of the estimation
histogram(result,20);    %the error of each iteration's estimation
xl=xlabel('Error $X - \tilde{X}$');
set(xl,'Interpreter','latex');
yl=ylabel('Number of error terms');
set(xl,'Interpreter','latex');
title('Histogram of error terms in the estimation');