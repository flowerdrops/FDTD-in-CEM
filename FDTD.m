close all
clear all
%FDTD计算方法，结果得到port1和port2两个时域数据
%传输线两边放挡板
fmax=10e9;  d=0.05; w=0.03; L=0.1;
%设置挡板条件
    k_x1=1/3;    k_x2=2/3;    k_z1=1/2;    k_v1=1/4;    k_v2=3/4;
    
c=3e8;T=1/(2*fmax);t0=3*T;lemda_min=c/fmax;dz=lemda_min/20;dx=dz/2;
niu=4*pi*1e-7;eku=1/(36*pi)*1e-9;%关于eku和niu的大小千万不可以搞错，不然会导致若干次迭代后非常严重的误差
dt=dx/(2*c);Nz=round(L/dz);Nx=round(d/dx);
% [nnz,nnx]=meshgrid(1:(Nz+1),1:Nx); 


x1=d*k_x1;x2=d*k_x2;z1=L*k_z1;v1=L*k_v1;v2=L*k_v2;

nx1=round(x1/dx);  nx2=round(x2/dx); nz1=round(z1/dz);  nv1=round(v1/dz); nv2=round(v2/dz);    nw=round(w/dz);

%设定数组初始值
Ex=zeros(Nx,Nz+1);Hy=zeros(Nx,Nz);Ez=zeros(Nx+1,Nz);Ex1=zeros(Nx,1);Ex2=zeros(Nx,1);

for tt=0:4500
%     tt
    t=tt*dt;
    F=exp(-(t-1*t0)^2./T^2);

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

    for i=2:Nx
        for j=1:Nz
            Ez(i,j)=Ez(i,j)+dt./eku.*(Hy(i,j)-Hy(i-1,j))./dx;
        end
    end
    Ez(1,1:Nz)=0;
    Ez(Nx+1,1:Nz)=0;
    Ez(nx1,nz1:nz1+nw)=0;
    Ez(nx2,nz1:nz1+nw)=0;
    
           mesh(Ex)
           axis([1,Nz,1,Nx,-2,2])
%            mesh(nnx,nnz,Ex)
           pause(0.01)

    Ex1=Ex(:,2);
    Ex2=Ex(:,Nz);
end











