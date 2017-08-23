function [X] =  createLH(pars,numsamples)

x_LB = pars(:,1)';
x_UB = pars(:,2)';

numpars = size(pars,1);

B_mult = ones(numsamples,1);

X_LB = B_mult*x_LB;
X_UB = B_mult*x_UB;

X_hat = lhsdesign(numsamples,numpars);

X = X_LB +(X_UB - X_LB).*X_hat;

end