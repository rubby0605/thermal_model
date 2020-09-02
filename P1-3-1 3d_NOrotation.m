%% RK try
% constant
G=6.67e-11;      % m^3/kg/s^2
K=1.38e-23;      % Boltzmann constant
M=2.306518e21; % kg
mole=6.02e23;  % moleculer 
distance=10.11595804; % [AU]
dv=1;
den_atm=zeros(ftd,ffd);   % Density of the atmosphere 
lifetime=(distance/0.386)^2*(2.83e6); %[t]
Radii=764.3e3;   % [m] P.C. Thomas et al,,2007
mmass=44e-3;  % Oxygen
ve=sqrt(2*G*M/Radii);          % escape speed [m/s]
%% Probability of sin
qq=1;
n_sin=0;
for i=0:0.01:pi 
    aa=sin(i);
    aa=aa*10000;
    p_sin(n_sin+1:n_sin+1+ceil(aa))=i;
    n_sin=n_sin+ceil(aa);
    qq=qq+1;
end
n_sin=n_sin+1;
%% Probability of cos
qq=1;
n_cos=0;
for i=-pi/2:0.01:pi/2
    aa=cos(i)*1e4;
    p_cos(n_cos+1:n_cos+1+ceil(aa))=i;
    n_cos=n_cos+ceil(aa);
    qq=qq+1;
end
n_cos=n_cos+1;
%% Main Program
   aaa=zeros(3,1000);
   ii=1;
   s=zeros(3,1);
   h=0.1;         % delta time
for n=1:1000
thi=p_cos(ceil(rand*n_cos));     % random theta
fi=p_sin(ceil(rand*n_sin));      % random phi
   td=ceil((thi+pi/2)/pi*ftd);
   fd=ceil(fi/2/pi*ffd);
   W=1;
   mm=0;
   local_T=plotT(td,fd);
   
while(W>=0.01&&local_T>=60)
   alt=0;
   td=ceil((thi+pi/2)/pi*ftd);
   fd=ceil(fi/2/pi*ffd);
   if(fd==0)
       fd=ffd;
   end
   if(td==0)
       td=ftd;
   end
   local_T=plotT(td,fd);
   flyd=rand*pi-pi/2;
   xi=rand*pi-pi/2;
   vr=0;
   for mk=1:12
   vr=vr+rand*12/11.9;
   end
   stdv2=2*K*local_T/mmass*mole;
   vr=(vr+6)/12*2*sqrt(stdv2);
   s(1)=Radii*cos(thi)*sin(fi); % x
   s(2)=Radii*cos(thi)*cos(fi); % y
   s(3)=Radii*sin(thi);         % z
   v(1)=vr*cos(-flyd+thi)*sin(xi);
   v(2)=vr*cos(-flyd+thi)*cos(xi);
   v(3)=vr*sin(-flyd+thi);
   oa=G*M*s(:)/(Radii^1.5);
   while(Ri>=Radii&&ii<=1000)
   for i=1:3
%       k1=h*(v(i)-G*M*s(i)/(Radii^1.5)*h);
%       k2=h*(v(i)-G*M*(s(i)+k1)/(Radii^1.5)*1.5*h);
%       k3=h*(v(i)-G*M*(s(i)+k2)/(Radii^1.5)*1.5*h);
%       k4=h*(v(i)-G*M*(s(i)+k3)/(Radii^1.5)*2*h);
       s(i)=s(i)+v(i)*h-G*M*s(i)/(Radii^1.5)*(2*h^2/3)-oa(i)*h^2/6;
       v(i)=v(i)-G*M*s(i)/(Radii^1.5)*h;
       oa(i)=G*M*s(i)/(Radii^1.5);
   end
   aaa(1,ii)=s(1);  %x
   aaa(2,ii)=s(2);  %y
   aaa(3,ii)=s(3);  %z
   Ri=sum(sqrt(s(:).^2));%%
   ii=ii+1;
   if(fi>=pi)
   alt=alt+h;
   W=1*exp(-alt/lifetime);
   end
   end
   thi=acos(s(3)/Ri);
   fi=atan(s(1)/s(2));
   if(thi<=-pi/2)
       if(fi>=pi/2)
       fi=fi-pi;
       thi=thi+pi/2;
       else
       fi=fi+2*pi;
       thi=-thi;
       end
   end
   if(thi>pi/2)
       if(fi>=pi/2)
       fi=-fi+pi;
       thi=-thi+pi;
       else
       fi=-fi+2*pi;
       thi=-thi+pi;
       end
   end
   if(fi<=0)
       fi=fi+2*pi;
   elseif(fi>2*pi)
       fi=fi-2*pi;
   end
end
end
%%

for i=1:1000
plot3(aaa(1,i),aaa(2,i),aaa(3,i))
hold on;
end