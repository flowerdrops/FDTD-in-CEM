%FDTD计算方法，结果得到port1和port2两个时域数据
%传输线两边放挡板
T=0.5e-9;
d=0.18; w=0.03;
L=1.2;
%设置挡板条件
    k_x1=1/3;
    k_x2=2/3;
    k_z1=1/2;
    k_v1=1/4;
    k_v2=3/4;
c=3e8;
t0=3*T;
fmax=1/(2*T);
lemda_min=c/fmax;
dz=lemda_min/40;
dx=dz/2;
niu=4*pi*1e-7;
eku=1/(36*pi)*1e-9;%关于eku和niu的大小千万不可以搞错，不然会导致若干次迭代后非常严重的误差
dt=dx/(2*c);
Nz=round(L/dz);
Nx=round(d/dx);
 
x1=d*k_x1;
x2=d*k_x2;
z1=L*k_z1;
v1=L*k_v1;
v2=L*k_v2;
 
nx1=round(x1/dx);  
nx2=round(x2/dx); 
nz1=round(z1/dz);
nv1=round(v1/dz);
nv2=round(v2/dz);
nw=round(w/dz);
%设定数组初始值
Ex=zeros(Nx,Nz+1);
Hy=zeros(Nx,Nz);
Ez=zeros(Nx+1,Nz);
Ex1=zeros(Nx,1);
Ex2=zeros(Nx,1);
port1=[];
port2=[];
for tt=0:4500
%     tt
    t=tt*dt;
    F=exp(-(t-2*t0)^2./T^2);
    
    
%     if tt<120
%         Ex(1:Nx,1)=F;
%     end
   %设定挡板条件，由于是金属，统统短路
    
    %中间放挡板
%     Ex(nx1:nx2,nz1)=0;
    %两边放挡板
%     Ex(1:nx1,nz1)=0;
%     Ex(nx2+1:Nx,nz1)=0;
%         Ex(nx1+1:nx2,nz1:nz1+nw)=0;
    %放满
%     Ex(1:Nx,nz1)=0;
    %设立观测点
    V1=0;%左边的
    V2=0;%右边的
    for i=1:Nx
        for j=1:Nz
            Hy(i,j)=Hy(i,j)+dt./niu.*((Ez(i+1,j)-Ez(i,j))./dx+(Ex(i,j)-Ex(i,j+1))./dz);
        end
    end
    
    for i=1:Nx
        for j=2:Nz
            if j==4
                Ex(i,j)=Ex(i,j)+dt./eku.*(Hy(i,j-1)-Hy(i,j))./dz+F;
            else
                Ex(i,j)=Ex(i,j)+dt./eku.*(Hy(i,j-1)-Hy(i,j))./dz;
            end
        end
    end
     %吸收边界
    Ex(1:Nx,1)=Ex1(1:Nx)+(c*dt-dz)/(c*dt+dz)*(Ex(1:Nx,2)-Ex(1:Nx,1));
    Ex(1:Nx,Nz+1)=Ex2(1:Nx)+(c*dt-dz)/(c*dt+dz)*(Ex(1:Nx,Nz)-Ex(1:Nx,Nz+1));
%     Ex(nx1+1:nx2,nz1:nz1+nw)=0;
%     Ex(nx1:nx2,nz1)=0;
%     Ex(nx1:nx2,nz1+nw)=0;
    for i=2:Nx
        for j=1:Nz
            Ez(i,j)=Ez(i,j)+dt./eku.*(Hy(i,j)-Hy(i-1,j))./dx;
        end
    end
    Ez(1,1:Nz)=0;
    Ez(Nx+1,1:Nz)=0;
    Ez(nx1,nz1:nz1+nw)=0;
    Ez(nx2,nz1:nz1+nw)=0;
%    Ez(nx1+1,nz1:nz1+nw)=0;
%    Ez(nx2+1,nz1:nz1+nw)=0;
%     Ez(nx1+1,nz1:nz1+nw)=0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Ez(nx2,nz1:nz1+nw)=0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %开始计算观测点的电压，由于线上可能不同，采用积分方式
 for i=1:Nx
     V1=V1+Ex(i,nv1);
     V2=V2+Ex(i,nv2);
 end
 %保存数据
 port1=[port1,Ex(:,nv1)];
 port2=[port2,Ex(:,nv2)];
%        if  mod(tt,1)==0
% %            plot(Ex(2,:))
% %            axis([1,Nz,-1,1])
            mesh(Ex)
            axis([1,Nz,1,Nx,-2,2])
           pause(0.01)
%        end
    Ex1=Ex(:,2);
    Ex2=Ex(:,Nz);
end
function [S11,S21,f,dt]=makeS(T,L,d,w,l_1,l_2,l_z,p1,p2)
%有挡板的情况
% T=0.2e-10;%脉冲宽度
% d=0.005;%传输线宽度
% L=0.02;%传输线长度
% w=0.0048;%金属板厚度
% p1=1/4;%端口1位置
% p2=3/4;%端口2位置
% 挡板位置
% l_1=2/10;
% l_2=9/10;
% l_z=1/2;
[port1,port2,dt]=my_FDTD(T,L,d,w,l_1,l_2,l_z,p1,p2);
%没挡板的情况
port0=my_FDTD0(T,L,d,p1,p2);
% figure
% plot([1:length(port0)]*dt,port0)
V_in=sum(port0);
V_re=sum(port1)-sum(port0);
V_tr=sum(port2);
% plot([1:length(V_re)]*dt,V_re,'r');
% title('反射波是红的，透射波是绿的')
% hold on
% plot([1:length(V_tr)]*dt,V_tr,'g');
 
[V_re_f,f]=fuliye15_25(V_re,dt,1e7,T);
[V_in_f,f]=fuliye15_25(V_in,dt,1e7,T);
[V_tr_f,f]=fuliye15_25(V_tr,dt,1e7,T);
% find(f==(1/(2*T)));
% [V_re_f,f]=myfft(V_re,dt);
% [V_in_f,f]=myfft(V_in,dt);
% [V_tr_f,f]=myfft(V_tr,dt);
 
 
% figure
% plot(f,abs(V_re_f),'r')
% title('频域内的反射波红色，透射波绿色')
% hold on
% plot(f,abs(V_tr_f),'g')
 
% S11=abs(V_re_f)./abs(V_in_f);
% S21=abs(V_tr_f)./abs(V_in_f);
S11=abs((V_re_f)./(V_in_f));
S21=abs((V_tr_f)./(V_in_f));
% figure
% 
% plot(f(1:round(length(f)*k)),S11(1:round(length(f)*k)),'r')
% title('S11参数红色，S21参数绿色')
% hold on
% plot(f(1:round(length(f)*k)),abs(S21(1:round(length(f)*k))),'g')
% 
% figure
% plot(f(1:round(length(f)*k)),abs(S11(1:round(length(f)*k)).^2+S21(1:round(length(f)*k)).^2))
% axis([f(1),f(end),0,2])
end