
%Preamble
clc
clear all
clear workspace
%Parameters - Fixed:
gamma=9.12;
p1=0.17;
p2=0.55*10^(-5);
p3=11.9*10^(-7); 
r=0.032;
d0=1*(1.032*10^5);
mu2=9.12;
k=10^9;
a=100;
m=1*2395/365;
mu1=21.05;
Runs=2000;
%%%%%%%%
%Random vectors:
Randomvec1= 0.5 + (2-0.5).*rand(Runs,1);
Randomvec2=0.5 + (2-0.5).*rand(Runs,1);
Randomvec3=0.5 + (2-0.5).*rand(Runs,1);
Randomvec4= 0.5 + (2-0.5).*rand(Runs,1);
Randomvec5=0.5 + (2-0.5).*rand(Runs,1);
Randomvec6=0.5 + (2-0.5).*rand(Runs,1);
%Parameters that vary (scaling a random vector by their model chosen value):
p2var = (0.55*10^(-5)*(Randomvec1))';
d0var = ((1.032*10^5)*Randomvec2)';
mu2var= (9.12*(Randomvec3))';
p1var= (0.17*(Randomvec4))';
mu1var=(21.05*(Randomvec5))';
avar=(100*(Randomvec6))';

rbottom=zeros(1,length(d0var));
rtop=zeros(1,length(d0var));
randomr=zeros(1,length(d0var));


m=zeros(1,length(d0var));
% Creating the range for parameter r and determining the dose for a
% random chosen value from this range 
for i = 1:length(p2var)
    
     p1=p1var(i);
     p2=p2var(i);
     d0=d0var(i);
     a=avar(i);
      mu2=mu2var(i);
      mu1=mu1var(i);

          rbottom(i)=p2*d0/mu2;
          rtop(i)=p2*d0/mu2+p1/mu1;
          randomr(i)=rbottom(i)+ (rtop(i)-rbottom(i))*rand(1,1);
          m(i)=mu1*mu1*a*(p2*d0-mu2*randomr(i))/(mu1*mu2*randomr(i)-p2*d0*mu1-mu2*p1);
end

%Scaling for MMC dose

mscaled=365*m;
mls=50;
Molw=334;
conv=10^(6);
mingrams=mscaled/(conv/(Molw*mls));
J=linspace(1,length(Randomvec1),length(Randomvec2));
semilogy(J,mingrams);
xlabel('Index number of the hypothetical patient') 
ylabel('Upper bound value [mg]') 
%hold on



hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = d0var;
Y =p2var;
Z1 = rbottom;
Z2 = rtop;
figure;scatter3(X,Y,Z1,'filled'); hold on;
scatter3(X,Y,Z2,'filled')
line([X; X], [Y; Y], [Z1; Z2], 'Color', 'k','LineWidth',1.5);%,
xlabel('d_{0}') 
ylabel('p_{2}') 
zlabel('Intervals of parameter r') 

rdist=zeros(1,length(d0var));
for i=1:length(rbottom)
    
rdist(i)=rtop(i)-rbottom(i);
end

