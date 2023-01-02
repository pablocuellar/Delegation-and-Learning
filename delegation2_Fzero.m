clc
clear all

r = 0.01;
D = 1;
lambda = 0.02;

%%%
delta = exp(-r*D);
phi = 1-exp(-lambda*D);
%%%

%delta = 0.99;
alpha = 0.4;
%p =	2/3;
p = 0.85;
%p_P = 3/4;
p_P = 9/10;
%lambda = 0.02;
%phi = 1-exp(-lambda);

h = alpha*p + (1-alpha)*(1-p);
h_P = alpha*p_P + (1-alpha)*(1-p_P);

bar_theta = alpha*p/h;
underline_theta = alpha*(1-p)/(1-h);

bar_theta_A_P = alpha*p/h - (1-alpha*p/h);
underline_theta_A_P = alpha*(1-p)/(1-h) - (1-alpha*(1-p)/(1-h));

bar_theta_P_P = max(alpha*p_P/h_P - (1-alpha*p_P/h_P),0);
underline_theta_P_P = alpha*(1-p_P)/(1-h_P) - (1-alpha*(1-p_P)/(1-h_P));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_del = ((1-delta)*(alpha-(1-alpha))+delta*phi*alpha)/(1-delta*(1-phi));
W_never = h_P*bar_theta_P_P;
W_bal = h*((1-delta)*bar_theta_A_P + delta*phi*alpha)/(1-h*delta*(1-phi)-(1-h)*delta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ite = 10000;

max_W =-999; 
max_F = -999;
max_tau = -999;
max_T = -999;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j =0:ite-1

T = j/1000;
F = 0;

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

W(j+1) = num/den;
tau_vec(j+1) = tau;
T_vec(j+1) = T;



if W(j+1)>max_W
max_W =  W(j+1);
max_F = F;
max_tau = tau;
max_T = T;
end


  
end



j/ite
end

W_del
W_never
W_bal 
max_W
%max_T
%max_tau

W_del_vec = ones(ite,1)*W_del';
W_never_vec = ones(ite,1)*W_never';
W_bal_vec = ones(ite,1)*W_bal';


figure(1)
plot(T_vec, W, '-', T_vec, W_del_vec, ':',   T_vec, W_never_vec,   '--', T_vec, W_bal_vec,'-.','LineWidth',4.0)
title('Non-Aligned Delegation - W vs T ' )
xlabel('T') 
ylabel('W') 
legend({'w','w_{FD}', 'w_{ND}', 'w_{BD}' },'Location','southeast ')
ylim([0.1 alpha-0.05])

