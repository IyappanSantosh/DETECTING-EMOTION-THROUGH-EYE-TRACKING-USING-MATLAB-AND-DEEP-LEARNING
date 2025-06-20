function [p] = ARMA(fdata)
n = length(fdata);
phi = -0.5;
theta1 = 0.9;
sigma = 1;
g = ((1+(theta1*theta1)+(2*phi*theta1)*sigma))/(1-(phi*phi));
p(1) = 1;
p(2) = ((1+(theta1*phi))*(theta1+phi)*sigma)/(g*(1-(phi*phi)));
for j = 3:n
    p(j) = phi*fdata(j-1);
end
end

