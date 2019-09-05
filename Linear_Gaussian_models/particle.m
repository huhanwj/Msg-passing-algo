function result=particle(k,x)
    result=struct();
    % Calculating the global function for X_k %
    x.W=x.W_forward{1,k}+x.W_backward{1,k};
    x.Wm=x.W_forward{1,k}*x.m_forward{1,k}+x.W_backward{1,k}*x.m_backward{1,k};
    x.V=inv(x.W);
    x.m=x.V*x.Wm;
    result.x1=x.m(1,1)+sqrt(x.V(1,1))*randn(1);
    result.x2=x.m(2,1)+sqrt(x.V(2,2))*randn(1);
end
    
    
   
    