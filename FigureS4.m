tic
%Parameters:
mu1=21.05;
gamma=9.12;
p1=0.12;
p2=0.5*10^(-5);
p3=15.9*10^(-7); 
r=0.02;
d0=(1.032*10^5);
mu2=9.12;
k=10^9;
a=100;
m=2395/365;
mu1=21.05;
%Initial conditions:
Initial=[m/mu1+10 (0+10)  (d0/mu2+10) ]; 
y0 = Initial;
opt = odeset('AbsTol',1e-9,'RelTol',1e-6);

tspan1 =[0,100];


    
[t,y] = ode45(@(t,y) odefcn(t,y,a,r,mu1,m,mu2,p1,p2,p3,d0,gamma,k), tspan1, y0, opt);
plot(t,y(:,1),'m',t,y(:,2),'k',t,y(:,3)/10^3,'b')
s=plot(t,y(:,1),'m',t,y(:,2),'k',t,y(:,3)/10^3,'b');
set(s,'linewidth',3);
legend('MMC (M)','Tumor cells (T)','Effector cells (E)')
xlabel('Time (Days)', 'FontSize',20, 'FontWeight', 'bold');
ylabel('# Cells or drug dose [\mu M]', 'FontSize',20, 'FontWeight', 'bold');
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
xlim([0 100])
ylim([0 12])
fontsize(16,"points")
fontweight='bold';



toc
%ODEs system:

function dydt =odefcn(t,y,a,r,mu1,m,mu2,p1,p2,p3,d0,gamma,k)
dydt = zeros(3,1);
M=y(1);
T=y(2);
E=y(3);
dydt = [ -mu1*M+m;
  -T*p1*M/(M+a)+r*T*(1-T/k)-T*(p2*E); 
  gamma*(p1*T*M/(M+a))+E*(-mu2)+d0-p3*E*T];

end
