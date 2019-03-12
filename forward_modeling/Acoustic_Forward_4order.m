%2D Cartesian coordinates Forward Acoustic wave in tunnel coded by Yuxiao Ren.
%Code reprogramed from SH forward wave.
%PML based on the paper ���Բ�����ģ����pml���ձ߽������ĸĽ������� et al., 2009��
%CFS-CPML based on the paper: a.����CFS-CPML�߽紦���LOVE�沨���޲��ģ��?
%   b.Complex frequency shifted convolution PML for FDTD modelling of elastic waves
clear;
close all;
clc;

%*********************************************************************** 
%% fundamental parameters ��������
%***********************************************************************

V_0=1000; % reference media velocity=1000 kg/m^3

freq=20;   % central frequency of source in Hz

time_window=0.6;    %�ܼ���ʱ��
dt=2.0*1.0e-4;
nt=round(time_window/dt);           %����ʱ���������?

%%********************************************************************
%% mesh parameters �������?
%********************************************************************

dx=single(2);				%������
dy=dx;

%********************************************************************
%% PML parameters
%********************************************************************

npml=single(30);
m=single(4);
kappa_max=single(1);
alpha_max=single(0);
bj=single(1);

%********************************************************************
%% ���ʲ��� and observation mode
%********************************************************************
 load velocity.mat;

nx=single(400);
ny=single(600);

param1=velocity./1000;
% param1=1.5.*ones(nx,ny);
fsx=[single(100)];
fsy=[single(300),single(450)]; %local 
recv_x=single(npml+1);
recv_y=single(npml+1):10:single(ny-npml);
re_num=length(recv_y)
fs_num=length(fsy);

%********************************************************************
%% Sources ����Դ 
%********************************************************************

% Ricker wavelet�׿��Ӳ�
source=zeros(1,nt);  
for i=1:nt
     A=1;t=i*dt-1/freq;
     source(i)=A*(1-2*(pi*freq*t)^2)*exp(-(pi*freq*t)^2);
end
% source=zeros(1,nt);   %Ricker wavelet�׿��Ӳ�
% for i=1:nt
%      A=1;t=i*dt-1/freq;
%      source(i)=A*(1-2*(pi*freq*t)^2)*exp(-(pi*freq*t)^2);
% end
% figure(2);
% plot(source);

%%********************************************************************
%% computing middle parameters �����м����?
%********************************************************************
tic;

disp('Computing middle paramters for ...');
% disp('Tyz'); %�൱������������λ����y�����һ�׵���?
V_r = zeros(nx,ny+1);
V_r(:,1:ny)=param1; V_r(:,ny+1)=V_r(:,ny);
V=V_r.*V_0; 

MP_Tyz_ya=ones(nx,ny+1); 

sig_max=bj.*(m+1)./(150.*pi.*dx).*sqrt(max(max(V)));
sigy=sig_max.*repmat(([npml-2:-1:1 1:npml-2]./npml).^m,[nx,1]); %pml������1
kappay=1+(kappa_max-1).*([npml-2:-1:1 1:npml-2]./npml).^m; %pml������1
kappay=repmat(kappay,[nx,1]);
alphay=alpha_max.*([npml-2:-1:1 1:npml-2]./npml); %=0 %pml������1
alphay=repmat(alphay,[nx,1]);
MP_Tyz_ya(:,[3:npml ny-npml+1:ny-2])=1./kappay;
MP_Tyz_yb=exp(-((sigy./kappay)+alphay));
MP_Tyz_yc=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Tyz_yb-1);
clear rho_r Vs_r Vc_r rho Vs Vc sig_max sigy kappay alphay

disp('Txz');%�൱������������λ����x�����һ�׵���?
V_r=zeros(nx+1,ny);
V_r(1:nx,:)=param1; V_r(nx+1,:)=V_r(nx,:);
V=V_r.*V_0; 

MP_Txz_xa=ones(nx+1,ny);

sig_max=bj.*(m+1)./(150.*pi.*dy).*sqrt(max(max(V)));
sigx=sig_max.*repmat(([npml-2:-1:1 1:npml-2]'./npml).^m,[1,ny]); %pml������1
kappax=1+(kappa_max-1).*([npml-2:-1:1 1:npml-2]'./npml).^m; %pml������1
kappax=repmat(kappax,[1,ny]);
alphax=alpha_max.*([npml-2:-1:1 1:npml-2]'./npml); %=0 %pml������1
alphax=repmat(alphax,[1,ny]);
MP_Txz_xa([3:npml nx-npml+1:nx-2],:)=1./kappax; 
MP_Txz_xb=exp(-((sigx./kappax)+alphax));
MP_Txz_xc=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Txz_xb-1);
clear rho_r Vs_r Vc_r rho Vs Vc sig_max sigx kappay alphay

disp('Vz'); %�൱������������λ�ƹ���ʱ���һ�׵��������������ٶ�?
V_r=param1; 
V=V_r.*V_0; 

MP_Vz=dt.*V.^2;
MP_Vz_xa=ones(size(V_r)); 
MP_Vz_ya=ones(size(V_r));

sig_max=bj.*(m+1)./(150.*pi.*dx).*sqrt(max(max(V)));
sigx=sig_max.*repmat((([npml-1:-1:1 1:npml-1]'-0.5)./npml).^m,[1,ny]);
kappax=1+(kappa_max-1).*(([npml-1:-1:1 1:npml-1]'-0.5)./npml).^m;
kappax=repmat(kappax,[1,ny]);
alphax=alpha_max.*(([npml-1:-1:1 1:npml-1]'-0.5)./npml); %=0
alphax=repmat(alphax,[1,ny]);
MP_Vz_xa([2:npml nx-npml+1:nx-1],:)=1./kappax;
MP_Vz_xb=exp(-((sigx./kappax)+alphax));
MP_Vz_xc=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Vz_xb-1);

sig_max=bj.*(m+1)./(150.*pi.*dy).*sqrt(max(max(V)));
sigy=sig_max.*repmat((([npml-1:-1:1 1:npml-1]-0.5)./npml).^m,[nx,1]);
kappay=1+(kappa_max-1).*(([npml-1:-1:1 1:npml-1]-0.5)./npml).^m;
kappay=repmat(kappay,[nx,1]);
alphay=alpha_max.*(([npml-1:-1:1 1:npml-1]-0.5)./npml); %=0
alphay=repmat(alphay,[nx,1]);
MP_Vz_ya(:,[2:npml ny-npml+1:ny-1])=1./kappay; 
MP_Vz_yb=exp(-((sigy./kappay)+alphay));
MP_Vz_yc=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Vz_yb-1);
clear rho_r Vs_r Vc_r rho sig Vs Vc  
clear sig_max sigx sigy alphay kappay
% clear param1 param2 param3 param4

%********************************************************************
%% ��ʼ��
%********************************************************************

Tyz=single(zeros(nx,ny+1));
UV_Tyz_y=single(zeros(nx,2*(npml-2)));
Txz=single(zeros(nx+1,ny));
UV_Txz_x=single(zeros(2*(npml-2),ny));
Vz=single(zeros(nx,ny));
UV_Vz_x=single(zeros(2*(npml-1),ny));
UV_Vz_y=single(zeros(nx,2*(npml-1)));

ob_data_Vz=zeros(length(recv_x),length(recv_y),nt);
% frwd_field=zeros(nx-2*npml,ny-2*npml,nt);

%********************************************************************
%% computing forward wavefiled
%********************************************************************
myavi=VideoWriter('2D_forward_Acoustic_waves_4_20.avi','Motion JPEG AVI');
myavi.FrameRate=5;
open(myavi);
c1=9/8;c2=-1/24;
for j=single(1):single(nt)
	
    Vz(fsx,fsy)=Vz(fsx,fsy) + source(j); %obmode_1
    %Vz(fsx(2),fsy(2))=Vz(fsx(2),fsy(2)) + source(j); % obmode_2
%     Vz(fsx,fsy)=Vz(fsx,fsy) + repmat(source(:,j)',length(fsx),1); % obmode_3
    
	Txz(3:nx-1,:)=Txz(3:nx-1,:) +dt.*MP_Txz_xa(3:nx-1,:).*(c1.*(Vz(3:nx-1,:)-Vz(2:nx-2,:))+c2.*(Vz(4:nx,:)-Vz(1:nx-3,:)))./dx;
	UV_Txz_x=MP_Txz_xb.*UV_Txz_x + MP_Txz_xc.*(c1.*(Vz([3:npml nx-npml+2:nx-1],:)-Vz([2:npml-1 nx-npml+1:nx-2],:))+c2.*(Vz([4:npml+1 nx-npml+3:nx],:)-Vz([1:npml-2 nx-npml:nx-3],:)))./dx;
	Txz([3:npml nx-npml+2:nx-1],:)=Txz([3:npml nx-npml+2:nx-1],:) + dt.*UV_Txz_x;
	 Txz([2 nx],:)=Txz([3 nx-1],:);
	Txz([1 nx+1],:)=Txz([2 nx],:);
		
		
	Tyz(:,3:ny-1)=Tyz(:,3:ny-1) + dt.*MP_Tyz_ya(:,3:ny-1).*(c1.*(Vz(:,3:ny-1)-Vz(:,2:ny-2))+c2.*(Vz(:,4:ny)-Vz(:,1:ny-3)))./dy;
	UV_Tyz_y=MP_Tyz_yb.*UV_Tyz_y + MP_Tyz_yc.*(c1.*(Vz(:,[3:npml ny-npml+2:ny-1])-Vz(:,[2:npml-1 ny-npml+1:ny-2]))+c2.*(Vz(:,[4:npml+1 ny-npml+3:ny])-Vz(:,[1:npml-2 ny-npml:ny-3])))./dy;
	Tyz(:,[3:npml ny-npml+2:ny-1])=Tyz(:,[3:npml ny-npml+2:ny-1]) + dt.*UV_Tyz_y;
     Tyz(:,[2 ny])=Tyz(:,[3 ny-1]);
    Tyz(:,[1 ny+1])=Tyz(:,[2 ny]);
    
	Vz(2:nx-1,2:ny-1)=Vz(2:nx-1,2:ny-1) + MP_Vz(2:nx-1,2:ny-1).*MP_Vz_xa(2:nx-1,2:ny-1).*(c1.*(Txz(3:nx,2:ny-1)-Txz(2:nx-1,2:ny-1))+c2.*(Txz(4:nx+1,2:ny-1)-Txz(1:nx-2,2:ny-1)))./dx + MP_Vz(2:nx-1,2:ny-1).*MP_Vz_ya(2:nx-1,2:ny-1).*(c1.*(Tyz(2:nx-1,3:ny)-Tyz(2:nx-1,2:ny-1))+c2.*(Tyz(2:nx-1,4:ny+1)-Tyz(2:nx-1,1:ny-2)))./dy;
	UV_Vz_x=MP_Vz_xb.*UV_Vz_x + MP_Vz_xc.*(c1.*(Txz([3:npml+1 nx-npml+2:nx],:)-Txz([2:npml nx-npml+1:nx-1],:))+c2.*(Txz([4:npml+2 nx-npml+3:nx+1],:)-Txz([1:npml-1 nx-npml:nx-2],:)))./dx;
	Vz([2:npml nx-npml+1:nx-1],:)=Vz([2:npml nx-npml+1:nx-1],:) + MP_Vz([2:npml nx-npml+1:nx-1],:).*UV_Vz_x;
	UV_Vz_y=MP_Vz_yb.*UV_Vz_y + MP_Vz_yc.*(c1.*(Tyz(:,[3:npml+1 ny-npml+2:ny])-Tyz(:,[2:npml ny-npml+1:ny-1]))+c2.*(Tyz(:,[4:npml+2 ny-npml+3:ny+1])-Tyz(:,[1:npml-1 ny-npml:ny-2])))./dy;
	Vz(:,[2:npml ny-npml+1:ny-1])=Vz(:,[2:npml ny-npml+1:ny-1]) + MP_Vz(:,[2:npml ny-npml+1:ny-1]).*UV_Vz_y;
    
    ob_data_Vz(:,:,j)=Vz(recv_x,recv_y);
%     frwd_field(:,:,j)=Vz(npml+1:nx-npml,1+npml:ny-npml);
    
	if mod(j,50)==0 %��ͼ�Ͷ���
        figure(3);
        set(gcf,'outerposition',get(0,'screensize')); %ȫ��
        imagesc(squeeze(Vz(npml+1:nx-npml,1+npml:ny-npml)));
%         line([1,tunnel_length-npml],[round((nx-tunnel_width)/2-npml),round((nx-tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % ����Ͻ���?
%         line([1,tunnel_length-npml],[round((nx+tunnel_width)/2-npml),round((nx+tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % ����½���?
%         line([tunnel_length-npml,tunnel_length-npml],[round((nx-tunnel_width)/2-npml),round((nx+tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % ������
%         line([301-npml,301-npml],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % ��ֱ����        
%         line([327,212],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 60����б����
%         line([392,252],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 45����б����
%         line([420,292],[73,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 30����б����
%         line([288,172],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 60����б����
%         line([342,212],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 45����б����
%         line([360,252],[52,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 30����б����
%         pos = [237 12 40 40]; rectangle('Position',pos,'Curvature',[1 1]);
        colorbar;
%          caxis([-1e-2,1e-2]);
        shading interp;
        grid on;
        axis equal tight
        set(gca,'fontsize',12,'fontweight', 'Bold');        
        h=suptitle(['Vz at step = ',num2str(j), ',  time = ',num2str(j*dt),' sec']);
        set(h,'fontsize',18,'fontweight', 'Bold');
         writeVideo(myavi,getframe(gcf));
	end
end
 close(myavi);
toc;

 save ob_data.mat ob_data_Vz
close all

time=0:dt:(length(ob_data_Vz)-1)*dt; %����ʱ����
r_num=size(ob_data_Vz,1);
r_num=re_num;
load('ob_data.mat');
datax1=squeeze(ob_data_Vz(1,:,:));
% datax2=squeeze(ob_data_Vz(2,:,:));
 datax11=datax1'; %datax21=datax2';
% AGCwindow=0.08;
% datax_AGC1=AGCgain(datax11,dt,AGCwindow,1); % AGC
% datax_AGC2=AGCgain(datax21,dt,AGCwindow,1); % AGC

h1=figure(1);
set(gcf,'outerposition',get(0,'screensize')); %ȫ��
subplot(1,2,1);
wigb(datax11,1,1:r_num,time)
% 	set(gcf,'Position',[200 50 543 650]); % �̶���ͼ�ߴ�
	set(gca,'FontName','Times New Roman','FontSize',18);
	title(['Original Seismic Record'],'FontName','Times New Roman','Fontsize',18,'Fontweight','Bold');
	xlabel('{{\it{Trace}}}','FontName','Times New Roman','FontSize',18);
	ylabel('{{\it{Time}}/s}','FontName','Times New Roman','FontSize',18);
