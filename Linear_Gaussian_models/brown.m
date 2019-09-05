% A simple Bronwinan Simulation %
function par = brown(k,step)
    rng('shuffle');
    par=struct();
    par.x1=cumsum(k*randn(step,1));
    par.x2=cumsum(k*randn(step,1));
   
end