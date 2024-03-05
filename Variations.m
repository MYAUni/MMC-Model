tic
%Parameters:
gamma=9.12;
p1=0.17;
p2=0.55*10^(-5);
p3=11.9*10^(-7); 
r=0.032;%Range:0.01-0.045
d0=1*(1.032*10^5);
mu2=9.12;
k=10^9;
a=100;
m=1*2395/365;
mu1=21.05;
%Initial conditions:
Initial=[ 0 1*10^7 d0/mu2 ]; %All scenarios: small- Initial(2)=5.3*10^6,
% Medium- Initial(2)=1*10^7 , Large Initial(2)=6.62*10^7 
y0 = Initial;
opt = odeset('AbsTol',1e-9,'RelTol',1e-6);


%With treatment:

tspan1 =[0,500];

d0=0.6*1.032*10^5;
 %r=0.01;   
 %p2=0.37*10^(-5);

[t,y] = ode45(@(t,y) odefcn(t,y,a,r,mu1,m,mu2,p1,p2,p3,d0,gamma,k), tspan1, y0, opt);

semilogy(t, y(:,2), 'Color', [0.4660 0.6740 0.1880]	); % Plotting with black color
hold on
s = semilogy(t, y(:,2), 'Color', [0.4660 0.6740 0.1880]	); % Plotting with black color
% Adjusting line width
set(s, 'LineWidth', 5);
% Define the y-axis tick locations
%ytick_values = 10.^8 .* (1:6); % [10^8, 2*10^8, 3*10^8, 4*10^8, 5*10^8, 6*10^8]

% Plot your data

% Set the y-axis ticks
%yticks(ytick_values);

% Adjust the y-axis tick labels
%yticklabels({'10^8', '2 \times 10^8', '3 \times 10^8', '4 \times 10^8', '5 \times 10^8', '6 \times 10^8'});

% Adjust other plot properties as needed
%grid on;

ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.YAxis.Exponent = 4;
xlim([0 400])
%ylim([0 4*10^8])
%ytick_values = 10.^8 .* (1:3); % [10^8, 2*10^8, 3*10^8, 4*10^8, 5*10^8, 6*10^8]
set(gca,'YminorTick','off')
%yticks([ 0 10^2 10^4 10^6 10^8 10^10 ])
%yticks([ 0 10 10^3 10^5 10^7 10^9 ])
xlabel('Time (Days)')
ylabel('# Cells')
fontsize(16,"points")
fontweight='bold';


hold on
%Without treatment:
y0 = Initial;
opt = odeset('AbsTol',1e-9,'RelTol',1e-6);

tspan1 =[0,500];

d0=1.032*10^5;
% r=0.025;   
%p2=0.45*10^(-5);

[t,y] = ode45(@(t,y) odefcn(t,y,a,r,mu1,m,mu2,p1,p2,p3,d0,gamma,k), tspan1, y0, opt);

semilogy(t, y(:,2), 'Color', [0.8500 0.3250 0.0980]	); % Plotting with black color
hold on
s = semilogy(t, y(:,2), 'Color', [0.8500 0.3250 0.0980]	); % Plotting with black color
% Adjusting line width
set(s, 'LineWidth', 5);



hold on
%Without treatment:
y0 = Initial;
opt = odeset('AbsTol',1e-9,'RelTol',1e-6);

tspan1 =[0,500];

d0=1.4*1.032*10^5;
% r=0.045;   
%p2=0.55*10^(-5);
[t,y] = ode45(@(t,y) odefcn(t,y,a,r,mu1,m,mu2,p1,p2,p3,d0,gamma,k), tspan1, y0, opt);

semilogy(t, y(:,2), 'Color', [0.4940 0.1840 0.5560]	); % Plotting with black color
hold on
s = semilogy(t, y(:,2), 'Color', [0.4940 0.1840 0.5560]		); % Plotting with black color
% Adjusting line width
set(s, 'LineWidth', 5);




hold off
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