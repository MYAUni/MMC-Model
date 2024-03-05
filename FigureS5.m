tic
%Parameters:
gamma=9.12;
p1=0.12;
p2=0.37*10^(-5);
p3=15.9*10^(-7); 
r=0.045;
d0=0.8*(1.032*10^5);
mu2=9.12;
k=10^9;
a=100;
m=1*2395/365;
mu1=21.05;
%Initial conditions:
Initial=[ 0 1*10^8 d0/mu2 ]; %For other scenarios: Initial(2)=4*10^6,  Initial(2)=1*10^7 
y0 = Initial;
opt = odeset('AbsTol',1e-9,'RelTol',1e-6);


%With treatment:

tspan1 =[0,500];


    
[t,y] = ode23s(@(t,y) odefcn(t,y,a,r,mu1,m,mu2,p1,p2,p3,d0,gamma,k), tspan1, y0, opt);

semilogy(t,y(:,2),'k',t,y(:,3),'b')
s=semilogy(t,y(:,2),'k',t,y(:,3),'b');
set(s,'linewidth',4);
set(s(2),'linewidth',4);
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.YAxis.Exponent = 4;
xlim([0 400])
ylim([0 3*10^9])
set(gca,'YminorTick','off')
%yticks([ 0 10^2 10^4 10^6 10^8 10^10 ])
%yticks([ 0 10 10^3 10^5 10^7 10^9 ])
legend('Tumor cells','Effector cells')
xlabel('Time (Days)')
ylabel('# Cells')
fontsize(16,"points")
fontweight='bold';


hold on
%Without treatment:
y0 = Initial;
opt = odeset('AbsTol',1e-9,'RelTol',1e-6);

tspan1 =[0,500];
m=10*2395/365;
mu1=21.05;
    
[t,y] = ode23s(@(t,y) odefcn(t,y,a,r,mu1,m,mu2,p1,p2,p3,d0,gamma,k), tspan1, y0, opt);
% Plot with specified line styles
s = semilogy(t, y(:,2), 'r', t, y(:,3), 'b');

% Set colors for each line
set(s(1), 'Color', [0.4940 0.1840 0.5560]	);  % Purple
set(s(2), 'Color', [0.4660 0.6740 0.1880]	);    % Dark Red

% Add legend
legend('Tumor cells', 'Effector cells');

% Set line widths
set(s(1), 'LineWidth', 5);  % Line 1
set(s(2), 'LineWidth', 6);  % Line 2

hold off
toc

%ODEs system:
function dydt =odefcn(t,y,a,r,mu1,m,mu2,p1,p2,p3,d0,gamma,k)
dydt = zeros(3,1);
M=y(1);
T=y(2);
E=y(3);
dydt = [ -mu1*M+m;
  -T*p1*M/(M+a)+r*T*(1-T/k)-T*(p2*E); %0.95*
  gamma*(p1*T*M/(M+a))+E*(-mu2)+d0-p3*E*T];%-p5*R*(D+E)/7
%*(1-E/10^9)
end
