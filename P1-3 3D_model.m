%% RK try

% constant
G=6.67e-11;      % m^3/kg/s^2
K=1.38e-23;      % Boltzmann constant
M=2.306518e21; % kg
mole=6.02e23;  % moleculer 
distance=10.11595804; % [AU]
dv=1;
lifetime=(distance/0.386)^2*(2.83e6); %[t]
Radii=764.3e3;   % [m] P.C. Thomas et al,,2007
mmass=44e-3;  % Oxygen
ve=sqrt(2*G*M/Radii);          % escape speed [m/s]
%% Read the Temperature file
fid=fopen('RHEA_T_NOVP.dat','r');
for ii=1:2
    line = fgetl(fid);
    %str=sscanf(line,'%d');
    str=sscanf(line, '%f');
    if(ii==1),ay=str;end
    if(ii==2),by=str;end
end
ftd=ay;
ffd=by;
plotT=zeros(ftd,ffd);
for ii=1:ftd
    for jj=1:ffd
    line = fgetl(fid);
    str=sscanf(line, '%f');
    plotT(ii,jj)=str;
    end
end
fclose(fid);
hold on;
imagesc(plotT)
colorbar

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

%% 
hd=100;   %% height grid points
maxheight=1e6;  %% [meter]
den_atm=zeros(hd,ftd,ffd);
h=10;    %% delta time
horizon_CO2=zeros(ftd,ffd);
%% Calculate for many times

for n=1:10000
   alt=0;
   thi=p_cos(ceil(rand*n_cos));     % random theta
   fi=p_sin(ceil(rand*n_sin));      % random phi
   td=ceil((thi+pi/2)/pi*ftd);
   fd=ceil(fi/2/pi*ffd);
   horizon_CO2(td,fd)=horizon_CO2(td,fd)-1;
   if(fd==0)
       fd=ffd;
   end
   if(td==0)
       td=ftd;
   end
   for step=1:1000
   local_T=plotT(td,fd);
   if(local_T<=60),horizon_CO2(td,fd)=horizon_CO2(td,fd)+1;break;end
   ii=1;
   flyd=rand*pi-pi/2;
   xi=rand*pi-pi/2;
   vr=0;
   for mk=1:12,vr=vr+rand*12/11.9;end
   stdv2=2*K*local_T/mmass*mole;
   vr=(vr)/12*2*sqrt(stdv2);
   s(1)=Radii*cos(thi)*sin(fi); % x
   s(2)=Radii*cos(thi)*cos(fi); % y
   s(3)=Radii*sin(thi);         % z
   v(1)=vr*cos(-flyd+thi)*sin(xi);
   v(2)=vr*cos(-flyd+thi)*cos(xi);
   v(3)=vr*sin(-flyd+thi);
   Ri=sqrt(sum(s(:).^2));%%
   oa=G*M*s(:)/(Ri^3);
   h=10;
   for ii=1:1000
   for i=1:3
       ns(i)=s(i)+v(i)*h-G*M*s(i)/(Ri^3)*(1.5*(h^2))-oa(i)*(h^2)/6;
       v(i)=v(i)-G*M*s(i)/(Ri^3)*h;
       oa(i)=G*M*s(i)/(Ri^3);
   end
   s(:)=ns(:);
   if(s(1)>=0),alt=alt+h;end
   W=exp(-alt/lifetime);
   if(W<=0.01),disp(W);break;end
   Ri=sqrt(sum(s(1:3).^2));
   if(Ri<=Radii),break;end
   thi=asin(s(3)/Ri);
   fi=acot(s(1)/s(2));
   if(s(1)<0),fi=fi+pi;end
   if(s(1)>0&&s(2)<0),fi=fi+2*pi;end
   td=ceil((thi+pi/2)/pi*ftd);
   fd=ceil(fi/2/pi*ffd);
   ht=ceil((Ri-Radii)/(maxheight-Radii)*h);
   if(fd==0),fd=ffd;end
   if(td==0),fd=ffd;end
   if(ht>100),disp('too height>"<');break;end
   den_atm(ht,td,fd)=den_atm(ht,td,fd)+W;
   end
   end
   if(ceil(n/100)-n/100==0),disp(n);end
end
%%
figure(3)
subplot(2,1,1);imagesc(plotT)
colorbar
%%
figure(1)
plot_den=zeros(hd,ftd);
plot_den(:,:)=den_atm(:,:,ceil(ffd/4));
imagesc(plot_den);
hold on;
colorbar
gcapoint1=[1 ceil(Radii/maxheight*100) 100];
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'ytick',gcapoint1);
set(gca,'xtick',gcapoint2);
set(gca,'xticklabel',{'-90','0','90'});
set(gca,'yticklabel',{'core','horizon','2e6[m]'});
xlabel('latitude');
ylabel('distance');
title('Local Time = 12:00');
%%
figure(2)
plot_den=zeros(hd,ffd);
plot_den(:,:)=den_atm(:,ceil(ftd/2),:);
imagesc(plot_den);
hold on;
colorbar
gcapoint1=[1 ceil(Radii/maxheight*100) 100];
gcapoint2=[1 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
set(gca,'ytick',gcapoint1);
set(gca,'xtick',gcapoint2);
set(gca,'xticklabel',{'6','12','18','24','6'});
set(gca,'yticklabel',{'core','horizon','2e6[m]'});
xlabel('Local time');
ylabel('distance');
title('thi=0')
%%
figure(3)
imagesc(horizon_CO2);
hold on;
colorbar
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'ytick',gcapoint2);
set(gca,'yticklabel',{'-90','0','90'});
gcapoint1=[1 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
set(gca,'xtick',gcapoint);
set(gca,'xticklabel',{'6','12','18','24','6'});
xlabel('Local Time');
ylabel('latitude');
title('horizontal CO2')
%%
subplot(2,1,2);
for td=1:ftd 
    for fd=1:ffd
        if(abs(plotT(td,fd)-60)<=1)
            hold on;
            plot(td,fd,'k');
        end
    end
end
%% 3-D model
figure(5)
    [x,y,z]=sphere;
    sphere(Radii.*x,Radii.*y,Radii.*z);
    hold on
for hh=1:hd
    Ri=(maxheight-Radii)/9*hh;
    for fd=0:ffd-1
        fi=fd/ffd*2*pi;
        for td=0:ftd-1
            thi=td/ftd*pi-pi/2;
            nn=ceil(den_atm(hh,td+1,fd+1));
            for n=1:nn
                fi=rand*dfi+fi;
                thi=rand*dthi+thi;
                xx=Ri*cos(thi)*cos(fi);
                yy=Ri*cos(thi)*sin(fi);
                zz=Ri*sin(thi);
                plot3(xx,yy,zz,'r')
                hold on
            end
        end
    end
end
