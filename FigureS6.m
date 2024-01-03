tic
%Parameters
mu1=21.05;
gamma=9.12;
p1=0.12;
p2=0.55*10^(-5);
p3=11.9*10^(-7); 
r=0.032;
d0=1*(1.032*10^5);
mu2=9.12;
k=10^9;
a=100;
%Initial conditions

%Case1
Initial=[ 0 2.15*10^7 1*10^3 ];
y0 = Initial;
opt = odeset('AbsTol',1e-9,'RelTol',1e-6);

tspan1 =[0,600];

m=2994/365;
    mu1=21.05; 


    
[t,y] = ode45(@(t,y) odefcn(t,y,a,r,mu1,m,mu2,p1,p2,p3,d0,gamma,k), tspan1, y0, opt);

semilogy(t,y(:,2))
s=semilogy(t,y(:,2));
set(s,'linewidth',5);
s(1).Color= [0.5 0 0.8];
xlim([0 600])
ylim([0 3*10^7])
set(gca,'YminorTick','off')
yticks([ 0  1*10^7 2*10^7  3*10^7 4*10^7])

legend('$T(0)=2.95\times10^{7}$')
xlabel('Time (Days)')
ylabel('# Cells')
fontsize(16,"points")
fontweight='bold';
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

hold on 

%Case 2

tspan1 =[0,600];
Initial=[ 0 2.14*10^7 1*10^3 ]; %Simulations for proposal: [ 0 3.7*10^7 1*10^7 ];
y0 = Initial;
opt = odeset('AbsTol',1e-9,'RelTol',1e-6);


    
[t,y] = ode45(@(t,y) odefcn(t,y,a,r,mu1,m,mu2,p1,p2,p3,d0,gamma,k), tspan1, y0, opt);

semilogy(t,y(:,2))
s=semilogy(t,y(:,2));
set(s,'linewidth',5);
s(1).Color= [0.6350 0.0780 0.1840];
legend('$T(0)=2.96\times10^{7}$')




%Case 3:
hold on 
tspan1 =[0,600];
Initial=[ 0 2.13*10^7 1*10^3 ]; 
y0 = Initial;
opt = odeset('AbsTol',1e-9,'RelTol',1e-6);

    
[t,y] = ode45(@(t,y) odefcn(t,y,a,r,mu1,m,mu2,p1,p2,p3,d0,gamma,k), tspan1, y0, opt);

semilogy(t,y(:,2))
s=semilogy(t,y(:,2));
set(s,'linewidth',5);
s(1).Color= [0 .7 .7];
legend('$T(0)=2.97\times10^{7}$')



hold off
toc

function dydt =odefcn(t,y,a,r,mu1,m,mu2,p1,p2,p3,d0,gamma,k)
dydt = zeros(3,1);
M=y(1);
T=y(2);
E=y(3);
dydt = [ -mu1*M+m;
  -T*p1*M/(M+a)+r*T*(1-T/k)-T*(p2*E); 
  gamma*(p1*T*M/(M+a))+E*(-mu2)+d0-p3*E*T];
end
