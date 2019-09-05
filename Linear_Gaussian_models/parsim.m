function output = parsim(k,x)
    output=struct();
    % Calculating the global function for X_k %
    x.W=x.W_forward{1,k}+x.W_backward{1,k};
    x.Wm=x.W_forward{1,k}*x.m_forward{1,k}+x.W_backward{1,k}*x.m_backward{1,k};
    x.V=inv(x.W);
    x.m=x.V*x.Wm;
    output.x1=x.m(1,1)+sqrt(x.V(1,1))*randn(1);
    output.vx1=x.m(2,1)+sqrt(x.V(2,2))*randn(1);
    output.x2=x.m(3,1)+sqrt(x.V(3,3))*randn(1);
    output.vx2=x.m(4,1)+sqrt(x.V(4,4))*randn(1);
end