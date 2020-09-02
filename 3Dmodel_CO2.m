clear all;
%% RK try
% constant
G=6.67e-11;      % m^3/kg/s^2
K=1.38e-23;      % Boltzmann constant
M=2.306518e21; % kg
mole=6.02e23;  % moleculer 
distance=10.11595804; % [AU]
dv=1;
lifetime=4.98e7; %[t]
Radii=764.3e3;   % [m] P.C. Thomas et al,,2007
mmass=44e-3;  % Oxygen
ve=sqrt(2*G*M/Radii);          % escape speed [m/s]
fprintf('Escape Velocity= %f\n',ve);
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
figure(10)
imagesc(plotT)
colorbar
dthi=pi/ftd;                   % thita的間隔
dfi=2*pi/ffd;                    % fi的間隔

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
den_atm_O=den_atm;
%% 
hd=50;   %% height grid points
maxheight=1.2e6;  %% [meter]
den_atm=zeros(hd,ftd,ffd);
h=50;    %% delta time
op=0;
horizon=zeros(ftd,ffd);
cc=0;
%% Calculate for many times
%figure(4)
aad=0;
bbd=0;
for n=1:100000
   W=1;
   thi=p_cos(ceil(rand*n_cos));     % random theta
   fi=p_sin(ceil(rand*n_sin));      % random phi
   td=ceil((thi+pi/2)/pi*ftd);
   fd=ceil(fi/2/pi*ffd);
%   td=ceil(rand*ftd);
%   fd=ceil(rand*ffd);
%   thi=td/ftd*pi;
%   fi=fd/ffd*2*pi;
   if(fd==0)
       fd=ffd;
   end
   if(td==0)
       td=ftd;
   end
   for step=1:1e7
   local_T=plotT(td,fd);
   if(local_T<=60),horizon(td,fd)=horizon(td,fd)+W;break;end
   ii=1;
   flyd=rand*pi/2;
   xi=rand*2*pi;
   if(abs(flyd)<=0.1),continue;end
   vr=0;
   for mk=1:12,vr=vr+rand*12/11.9;end
   stdv2=2*K*local_T/mmass*mole;
   vr=(vr)/12*2*sqrt(stdv2);
   vr=abs(vr);
   s(1)=Radii*cos(thi)*sin(fi); % x
   s(2)=Radii*cos(thi)*cos(fi); % y
   s(3)=Radii*sin(thi);         % z
   v(1)=-sin(flyd)*sin(-xi)*cos(fi);
   y1=-sin(flyd)*cos(-xi);
   z1=cos(flyd);
   v(2)=(y1*sin(thi)+z1*cos(thi))*cos(fi);
   v(3)=(-y1*cos(thi)+z1*sin(thi));
   v(:)=v(:).*vr;
   Ri=sqrt(sum(s(:).^2));%%
   oa=G*M*s(:)/(Ri^3);
   for ii=1:1000
   for i=1:3
       ns(i)=s(i)+v(i)*h-G*M*s(i)/(Ri^3)*(1.5*(h^2))-oa(i)*(h^2)/6;
       v(i)=v(i)-G*M*s(i)/(Ri^3)*h;
       oa(i)=G*M*s(i)/(Ri^3);
   end
   s(:)=ns(:);
   if(s(1)>=0),W=W*exp(-h/lifetime);end
   if(W<=0.1),break;end
   Ri=sqrt(sum(s(1:3).^2));
   if(Ri<=Radii),break;end
   thi=asin(s(3)/Ri);
   fi=asin(s(1)/sqrt((s(1)^2+s(2)^2)));
   if(abs(s(1))<=0.1&&abs(s(2))<=0.1),fi=0;end
   if(s(1)>0&&s(2)<0), fi=-fi+pi;end
   if(s(1)<0&&s(2)<0), fi=-fi+pi;end
   if(s(1)<0&&s(2)>0), fi=fi+2*pi;end
   if(abs(thi)-abs(thi)~=0||abs(thi)<0),disp('??');break;end
   td=ceil((thi+pi/2)/pi*ftd);
   fd=ceil(fi/2/pi*ffd);
   ht=ceil((Ri-Radii)/(maxheight-Radii)*hd);
   if(ht<=0),ht=1;end
   if(fd==0),fd=ffd;end
   if(td==0),td=ftd;end
   if(ht>hd),disp('too height>"<');disp(ht);break;end
   den_atm(ht,td,fd)=den_atm(ht,td,fd)+W;
   if(W<=0.1),disp(W);break;end;
   end
   disp(n);
   
   if(bbd==fd),bbd=0;break;else bbd=fd;end;
   if(aad==td),aad=0;break;else aad=td;end;
   end
   if(bbd==0),continue;end;
   if(aad==0),continue;end;
end
cc=cc+n-1;

%% Normalize
toP=zeros(ftd,ffd);
toZ=zeros(ftd,ffd);
toP(:,:)=133.32237.*10.^(-1275.6./plotT(:,:)+0.00683*plotT(:,:)+8.307);
toZ(:,:)=sqrt(1./mmass./(2.*pi.*K.*plotT(:,:))).*toP(:,:).*1e4;  % [mole/cm^2/s];

for td=1:ftd
    for fd=1:ffd
        den_atm(:,td,fd)=den_atm(:,td,fd).*toZ(td,fd);
    end
end


%%
figure(2)

plot_den=zeros(hd,ftd);
plot_den(:,:)=den_atm(:,:,ceil(ffd/2));
imagesc(plot_den);
hold on;
colorbar
gcapoint1=[1 hd/2 hd];
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'ytick',gcapoint1);
set(gca,'xtick',gcapoint2);
set(gca,'xticklabel',{'-90','0','90'});
RR=(maxheight-Radii)/1000;
RR2=RR/2;
str=num2str(RR);
str2=num2str(RR2);
set(gca,'yticklabel',{'horizon',str2,str});
xlabel('latitude');
ylabel('distance [km]');
title('Local Time = 12:00');
%%
figure(3)
plot_den=zeros(hd,ffd);
plot_den(:,:)=den_atm(:,ceil(ftd/4),:);
imagesc(plot_den);
hold on;
colorbar
gcapoint1=[1 hd/2 hd];
gcapoint2=[1 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
set(gca,'ytick',gcapoint1);
set(gca,'xtick',gcapoint2);
set(gca,'xticklabel',{'6','12','18','24','6'});
set(gca,'yticklabel',{'horizon',str2,str});
xlabel('Local time');
ylabel('height [km]');
title('thi=0.3977')
%%
figure(4)
imagesc(horizon);
hold on;
colorbar
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'ytick',gcapoint2);
set(gca,'yticklabel',{'-90','0','90'});
gcapoint1=[1 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
set(gca,'xtick',gcapoint1);
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
%% 2-D model
figure(5)
        fd=ceil(ffd/4);
        fi=fd/ffd*2*pi;
for td=0:ftd-1
            thi=td/ftd*pi-pi/2;
            xx=Ri*cos(thi)*cos(fi);
            zz=Ri*sin(thi);
            plot(xx,zz,'b')
            hold on
end
        fd=ceil(ffd/4*3);
        fi=fd./ffd.*2.*pi;
for td=0:ftd-1
            thi=td/ftd*pi-pi/2;
            xx=Ri*cos(thi)*cos(fi);
            zz=Ri*sin(thi);
            plot(xx,zz,'b')
            hold on
end
for hh=1:hd
    Ri=(maxheight-Radii)/hd*hh+Radii;
            fd=ceil(ffd/4);
        fi=fd/ffd*2*pi;
    for td=0:ftd-1
            thi=td/ftd*pi-pi/2;
            if(den_atm(hh,td+1,fd+1)<=5),continue;end
            nn=ceil(den_atm(hh,td+1,fd+1));
            for n=1:nn
                fi=rand*dfi+fi;
                thi=rand*dthi+thi;
                xx=Ri*cos(thi);
                yy=Ri*cos(thi);
                zz=Ri*sin(thi);
                plot(xx,zz,'r')
                hold on
            end
    end
    disp(hh)
end
for hh=1:hd
    Ri=(maxheight-Radii)/hd*hh+Radii;
            fd=ceil(ffd/4);
        fi=fd/ffd*2*pi;
    for td=0:ftd-1
            thi=td/ftd*pi-pi/2;
            if(den_atm(hh,td+1,fd+1)<=5),continue;end
            nn=ceil(den_atm(hh,td+1,fd+1));
            for n=1:nn
                fi=rand*dfi+fi;
                thi=rand*dthi+thi;
                xx=Ri*cos(thi);
                yy=Ri*cos(thi);
                zz=Ri*sin(thi);
                plot(xx,zz,'r')
                hold on
            end
    end
    disp(hh)
end

    %%
    hold on
for hh=1:hd
    Ri=(maxheight-Radii)/hd*hh+Radii;
    for fd=0:ffd-1
        fi=fd/ffd*2*pi;
        for td=0:ftd-1
            thi=td/ftd*pi-pi/2;
            if(den_atm(hh,td+1,fd+1)<=5),continue;end
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
    disp(hh)
end
