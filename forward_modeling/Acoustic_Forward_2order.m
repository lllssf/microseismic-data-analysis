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

time_window=30;    %�ܼ���ʱ��
dt=1.0*1.0e-3;
nt=round(time_window/dt);           %����ʱ���������?

%%********************************************************************
%% mesh parameters �������?
%********************************************************************

dx=single(10);				%������
dy=dx;

%********************************************************************
%% PML parameters
%********************************************************************

npml=single(20);
m=single(4);
kappa_max=single(1);
alpha_max=single(0);
bj=single(1);

%********************************************************************
%% ���ʲ��� and observation mode
%********************************************************************
load velocity.mat;

nx=single(360);
ny=single(1280);

param1=velocity./1000;
fsx=single(npml+1);
fsy=single(640);
recv_x=single(npml+1);
recv_y=single(npml+1):10:single(ny-npml);
re_num=length(recv_y);
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
% figure(2);
% plot(source);

%%********************************************************************
%% computing middle parameters �����м����?
%********************************************************************
tic;

disp('Computing middle paramters for ...');
disp('Tyz'); 
V_r = zeros(nx,ny+1);
V_r(:,1:ny)=param1; V_r(:,ny+1)=V_r(:,ny);
V=V_r.*V_0; 

MP_Tyz_ya=ones(nx,ny+1); 

sig_max=bj.*(m+1)./(150.*pi.*dx).*sqrt(max(max(V)));
sigy=sig_max.*repmat(([npml-1:-1:1 1:npml-1]./npml).^m,[nx,1]); %pml������1
kappay=1+(kappa_max-1).*([npml-1:-1:1 1:npml-1]./npml).^m; %pml������1
kappay=repmat(kappay,[nx,1]);
alphay=alpha_max.*([npml-1:-1:1 1:npml-1]./npml); %=0 %pml������1
alphay=repmat(alphay,[nx,1]);
MP_Tyz_ya(:,[2:npml ny-npml+2:ny])=1./kappay;
MP_Tyz_yb=exp(-((sigy./kappay)+alphay));
MP_Tyz_yc=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Tyz_yb-1);
clear rho_r Vs_r Vc_r rho Vs Vc sig_max sigy kappay alphay

disp('Txz');%�൱������������λ����x�����һ�׵���?
V_r=zeros(nx+1,ny);
V_r(1:nx,:)=param1; V_r(nx+1,:)=V_r(nx,:);
V=V_r.*V_0; 

MP_Txz_xa=ones(nx+1,ny);

sig_max=bj.*(m+1)./(150.*pi.*dy).*sqrt(max(max(V)));
sigx=sig_max.*repmat(([npml-1:-1:1 1:npml-1]'./npml).^m,[1,ny]); %pml������1
kappax=1+(kappa_max-1).*([npml-1:-1:1 1:npml-1]'./npml).^m; %pml������1
kappax=repmat(kappax,[1,ny]);
alphax=alpha_max.*([npml-1:-1:1 1:npml-1]'./npml); %=0 %pml������1
alphax=repmat(alphax,[1,ny]);
MP_Txz_xa([2:npml nx-npml+2:nx],:)=1./kappax; 
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
sigx=sig_max.*repmat((([npml:-1:1 1:npml]'-0.5)./npml).^m,[1,ny]);
kappax=1+(kappa_max-1).*(([npml:-1:1 1:npml]'-0.5)./npml).^m;
kappax=repmat(kappax,[1,ny]);
alphax=alpha_max.*(([npml:-1:1 1:npml]'-0.5)./npml); %=0
alphax=repmat(alphax,[1,ny]);
MP_Vz_xa([1:npml nx-npml+1:nx],:)=1./kappax;
MP_Vz_xb=exp(-((sigx./kappax)+alphax));
MP_Vz_xc=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Vz_xb-1);

sig_max=bj.*(m+1)./(150.*pi.*dy).*sqrt(max(max(V)));
sigy=sig_max.*repmat((([npml:-1:1 1:npml]-0.5)./npml).^m,[nx,1]);
kappay=1+(kappa_max-1).*(([npml:-1:1 1:npml]-0.5)./npml).^m;
kappay=repmat(kappay,[nx,1]);
alphay=alpha_max.*(([npml:-1:1 1:npml]-0.5)./npml); %=0
alphay=repmat(alphay,[nx,1]);
MP_Vz_ya(:,[1:npml ny-npml+1:ny])=1./kappay; 
MP_Vz_yb=exp(-((sigy./kappay)+alphay));
MP_Vz_yc=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Vz_yb-1);
clear rho_r Vs_r Vc_r rho sig Vs Vc  
clear sig_max sigx sigy alphay kappay
clear param1 param2 param3 param4

%********************************************************************
%% ��ʼ��
%********************************************************************

Tyz=single(zeros(nx,ny+1));
UV_Tyz_y=single(zeros(nx,2*(npml-1)));
Txz=single(zeros(nx+1,ny));
UV_Txz_x=single(zeros(2*(npml-1),ny));
Vz=single(zeros(nx,ny));
UV_Vz_x=single(zeros(2*npml,ny));
UV_Vz_y=single(zeros(nx,2*npml));

ob_data_Vz=zeros(length(recv_x),length(recv_y),nt);
% frwd_field=zeros(nx-2*npml,ny-2*npml,nt);

%********************************************************************
%% computing forward wavefiled
%********************************************************************
myavi=VideoWriter('2D_forward_Acoustic_waves.avi','Motion JPEG AVI');
myavi.FrameRate=5;
open(myavi);

for j=single(1):single(nt)
	
%     Vz(fsx,fsy)=Vz(fsx,fsy) + source(j); %obmode_1
%     Vz(fsx,fsy)=Vz(fsx,fsy) + source(:,j); % obmode_2
    Vz(fsx,fsy)=Vz(fsx,fsy) + repmat(source(:,j)',length(fsx),1); % obmode_3
    
	Txz(2:nx,:)=Txz(2:nx,:) + dt.*MP_Txz_xa(2:nx,:).*(Vz(2:nx,:)-Vz(1:nx-1,:))./dx;
	UV_Txz_x=MP_Txz_xb.*UV_Txz_x + MP_Txz_xc.*(Vz([2:npml nx-npml+2:nx],:)-Vz([1:npml-1 nx-npml+1:nx-1],:))./dx;
	Txz([2:npml nx-npml+2:nx],:)=Txz([2:npml nx-npml+2:nx],:) + dt.*UV_Txz_x;
	Txz([1 nx+1],:)=Txz([2 nx],:);
		
	Tyz(:,2:ny)=Tyz(:,2:ny) + dt.*MP_Tyz_ya(:,2:ny).*(Vz(:,2:ny)-Vz(:,1:ny-1))./dy;
	UV_Tyz_y=MP_Tyz_yb.*UV_Tyz_y + MP_Tyz_yc.*(Vz(:,[2:npml ny-npml+2:ny])-Vz(:,[1:npml-1 ny-npml+1:ny-1]))./dy;
	Tyz(:,[2:npml ny-npml+2:ny])=Tyz(:,[2:npml ny-npml+2:ny]) + dt.*UV_Tyz_y;
    Tyz(:,[1 ny+1])=Tyz(:,[2 ny]);
    
	Vz=Vz + MP_Vz.*MP_Vz_xa.*(Txz(2:nx+1,:)-Txz(1:nx,:))./dx + MP_Vz.*MP_Vz_ya.*(Tyz(:,2:ny+1)-Tyz(:,1:ny))./dy;	
	UV_Vz_x=MP_Vz_xb.*UV_Vz_x + MP_Vz_xc.*(Txz([2:npml+1 nx-npml+2:nx+1],:)-Txz([1:npml nx-npml+1:nx],:))./dx;
	Vz([1:npml nx-npml+1:nx],:)=Vz([1:npml nx-npml+1:nx],:) + MP_Vz([1:npml nx-npml+1:nx],:).*UV_Vz_x;
	UV_Vz_y=MP_Vz_yb.*UV_Vz_y + MP_Vz_yc.*(Tyz(:,[2:npml+1 ny-npml+2:ny+1])-Tyz(:,[1:npml ny-npml+1:ny]))./dy;
	Vz(:,[1:npml ny-npml+1:ny])=Vz(:,[1:npml ny-npml+1:ny]) + MP_Vz(:,[1:npml ny-npml+1:ny]).*UV_Vz_y;
    
%     ob_data_Vz(:,:,j)=Vz(recv_x,recv_y);
%     frwd_field(:,:,j)=Vz(npml+1:nx-npml,1+npml:ny-npml);
    
	if mod(j,50)==0 %��ͼ�Ͷ���
        figure(3);
        set(gcf,'outerposition',get(0,'screensize')); %ȫ��
        imagesc(squeeze(Vz(npml+1:nx-npml,1+npml:ny-npml)));
%         line([1,tunnel_length-npml],[round((nx-tunnel_width)/2-npml),round((nx-tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % ����Ͻ���?
%         line([1,tunnel_length-npml],[round((nx+tunnel_width)/2-npml),round((nx+tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % ����½���?
%         line([tunnel_length-npml,tunnel_length-npml],[round((nx-tunnel_width)/2-npml),round((nx+tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % ������
%         line([301-npml,301-npml],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % ��ֱ����        
%         line([400,165],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 60����б����
%         line([360,200],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 45����б����
%         line([360,141],[35,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 30����б����
%         pos = [237 12 40 40]; rectangle('Position',pos,'Curvature',[1 1]);
        colormap(gray);
        colorbar;
        shading interp;
        grid on;
        axis equal tight;
        set(gca,'fontsize',12,'fontweight', 'Bold');        
        h=suptitle(['Vz at step = ',num2str(j), ',  time = ',num2str(j*dt),' sec']);
        set(h,'fontsize',18,'fontweight', 'Bold');
        writeVideo(myavi,getframe(gcf));
	end
end
close(myavi);
toc;

% save ob_data.mat ob_data_Vz
close all

