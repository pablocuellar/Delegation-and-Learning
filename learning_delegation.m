clc
clear all

delta = 0.99;
alpha = 0.3;
p =	0.51;
p_P = 0.6;
lambda = 2;
phi = 1-exp(-lambda );

h = alpha*p + (1-alpha)*(1-p);
h_P = alpha*p_P + (1-alpha)*(1-p_P);

bar_theta = alpha*p/h;
underline_theta = alpha*(1-p)/(1-h);

bar_theta_A_P = alpha*p/h - (1-alpha*p/h);
underline_theta_A_P = alpha*(1-p)/(1-h) - (1-alpha*(1-p)/(1-h));

bar_theta_P_P = max(alpha*p_P/h_P - (1-alpha*p_P/h_P),0);
underline_theta_P_P = alpha*(1-p_P)/(1-h_P) - (1-alpha*(1-p_P)/(1-h_P));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_del = ((1-delta)*(alpha-(1-alpha))+delta*phi*alpha)/(1-delta*(1-phi))
W_never = h_P*bar_theta_P_P
W_bal = h*((1-delta)*bar_theta_A_P + delta*phi*alpha)/(1-h*delta*(1-phi)-(1-h)*delta)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ite = 10000;

max_W =0; 
max_F = -999;
max_tau = -999;
max_T = -999;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j =0:ite-1

T = j/1000;


max_WW = -9999; 
for i = 0:ite-1
    
F = i/1000;

V_star = (delta^F)*alpha;
x = delta*(1-phi);
V_N = (1/(1-delta*x^T))*  ( ((1-x^T)/(1-x)) * ((1-delta)*alpha + delta*phi*V_star) + x^T *h*(1-delta)*(bar_theta-underline_theta));
V = h*(1-delta)*(bar_theta - underline_theta) + delta*V_N;
V_I = V_N-((1-delta)/delta)*underline_theta;
tau = log((V_I-phi*V_star)/((1-phi)*V))   /  log(delta);


if tau >=0


W_star = (1-delta^F)*h_P*bar_theta_P_P + (delta^F)*alpha;
A_I = phi*W_star  + (1-phi)*((1-delta^tau)*h_P*bar_theta_P_P);
A_N = ((1-(delta*(1-phi))^T)/(1-(delta*(1-phi)))) * ((1-delta)*(alpha-(1-alpha)) + delta*phi*W_star) ;

num = h*(1-delta)*bar_theta_A_P + h*delta*A_I + (1-h)*delta*A_N;
den = 1-(delta^(tau+1))*h*(1-phi) - (delta^(T+1))*(1-h)*(1-phi)^T;

W(i+1) = num/den;

F_vec(i+1) = F;
tau_vec(i+1) = tau;


if W(i+1)>max_WW
max_WW =  W(i+1);
end


if W(i+1)>max_W
max_W =  W(i+1);
max_F = F;
max_tau = tau;
max_T = T;
end
    

end

%i/ite    
end

T_vec(j+1) = T;
opt_W(j+1) = max_WW;

%j/ite
end



max_W
max_T
max_F
max_tau



figure(1)
plot(T_vec, opt_W )

