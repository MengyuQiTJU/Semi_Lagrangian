clc;
clear;
tic;
X=12800;%//length of the domain
Y=X;
T=12000;%//time of simulation
Nx=128;%//the number of cells in space
Ny=Nx;
M=Nx/128*1200;%//the number of cells in time,time step
w=pi/6000;
dx=X/Nx;% //this is the size of the space step
dy=Y/Ny;
% sx=141.43*5;
sx=141.43*5;
sy=100*5;
dt=T/M; %//this is the size of the time step
D=5;%%%扩散系数
x1=zeros(1,Nx+1); % array allocation
y1=zeros(1,Ny+1);
c1=zeros(Nx+1,Ny+1);
cx1=zeros(Nx+1,Ny+1);
cy1=zeros(Nx+1,Ny+1);
cxy1=zeros(Nx+1,Ny+1);
u=zeros(Nx+1,Ny+1);
v=zeros(Nx+1,Ny+1);
r=zeros(Nx+1,Ny+1);
Cr_x=zeros(Nx+1,Ny+1);
Cr_y=zeros(Nx+1,Ny+1);
a_x=zeros(Nx+1,Ny+1);
a_y=zeros(Nx+1,Ny+1);
a1_x=zeros(Nx+1,Ny+1);
a2_x=zeros(Nx+1,Ny+1);
a3_x=zeros(Nx+1,Ny+1);
a4_x=zeros(Nx+1,Ny+1);
a1_y=zeros(Nx+1,Ny+1);
a2_y=zeros(Nx+1,Ny+1);
a3_y=zeros(Nx+1,Ny+1);
a4_y=zeros(Nx+1,Ny+1);
b1_x=zeros(Nx+1,Ny+1);
b2_x=zeros(Nx+1,Ny+1);
b3_x=zeros(Nx+1,Ny+1);
b4_x=zeros(Nx+1,Ny+1);
b1_y=zeros(Nx+1,Ny+1);
b2_y=zeros(Nx+1,Ny+1);
b3_y=zeros(Nx+1,Ny+1);
b4_y=zeros(Nx+1,Ny+1);
R=zeros(Nx+1,Ny+1);
seita=zeros(Nx+1,Ny+1);
xd=zeros(Nx+1,Ny+1);
yd=zeros(Nx+1,Ny+1);
show=zeros(Nx+1,Ny+1);
x0=0;
y0=3200;
rc=0.1;
m=1;
for i=1:Nx+1
    for j=1:Ny+1
    x1(i)=dx*(i-(Nx/2+1));
    y1(j)=dy*(j-(Ny/2+1));
    r(i,j)=sqrt((x1(i)-x0)^2+(y1(j)-y0)^2);
    c1(i,j)=exp(-(x1(i)-x0)^2/(2*sx^2)-(y1(j)-y0)^2/(2*sy^2));
    cx1(i,j)=-c1(i,j)*(x1(i)-x0)/sx^2;
    cy1(i,j)=-c1(i,j)*(y1(j)-y0)/sy^2;
%     if r(i,j)<=rc
%         c1(i,j)=5*(1+cos(r(i,j)*pi/rc));
%         cx1(i,j)=-(5*pi*sin((pi*r(i,j)/rc))*(x1(i)-x0))/(rc*r(i,j));
%         cy1(i,j)=-(5*pi*sin((pi*r(i,j)/rc))*(y1(j)-y0))/(rc*r(i,j));
%     else
%         c1(i,j)=0;
%         cx1(i,j)=0;
%         cy1(i,j)=0;
%     end

%     if isnan(cx1(i,j))
%         cx1(i,j)=0;
%     end
%     if isnan(cy1(i,j))
%         cy1(i,j)=0;
%     end

    u(i,j)=-w*y1(j);
    v(i,j)=w*x1(i);
    R(i,j)=sqrt(x1(i)^2+y1(j)^2);
    if i>Nx/2+1 % for quadrant I & IV & X positive half shaft
        seita(i,j)=atan(y1(j)/x1(i));
    elseif i<Nx/2+1 % for quadrant II & III & X negative half shaft
        seita(i,j)=atan(y1(j)/x1(i))+pi;
    elseif i==Nx/2+1 && j>Ny/2+1 % for Y positive half shaft
        seita(i,j)=pi/2;
    elseif i==Nx/2+1 && j<Ny/2+1 % for Y negative half shaft
        seita(i,j)=-pi/2;
    else
        seita(i,j)=0;
    end
    xd(i,j)=R(i,j)*cos(seita(i,j)-w*dt);
    yd(i,j)=R(i,j)*sin(seita(i,j)-w*dt);
    end
end
% cxt=zeros(Nx+1,Ny+1);

for i=2:Nx
    for j=2:Ny
        cxy1(i,j)=c1(i,j)*(x1(i)-x0)/sx^2*(y1(j)-y0)/sy^2;
    end
end

for i=2:Nx
    for j=2:Ny
        %     Cr_x(i,j)=abs(x1(i)-xd(i,j))/dx;
        %     Cr_y(i,j)=abs(y1(j)-yd(i,j))/dy;
        if i>=Nx/2+1 && j>Ny/2+1 % for quadrant I & Y positive half shaft
            Cr_x(i,j)=(x1(i+1)-xd(i,j))/dx;
            Cr_y(i,j)=(y1(j)-yd(i,j))/dy;
        elseif i<Nx/2+1 && j>=Ny/2+1 % for quadrant II & X negative half shaft
            Cr_x(i,j)=(x1(i+1)-xd(i,j))/dx;
            Cr_y(i,j)=(y1(j+1)-yd(i,j))/dy;
        elseif i<=Nx/2+1 && j<Ny/2+1 % for quadrant III & Y negative half shaft
            Cr_x(i,j)=(x1(i)-xd(i,j))/dx;
            Cr_y(i,j)=(y1(j+1)-yd(i,j))/dy;
        elseif i>Nx/2+1 && j<=Ny/2+1 % for quadrant IV & X positive half shaft
            Cr_x(i,j)=(x1(i)-xd(i,j))/dx;
            Cr_y(i,j)=(y1(j)-yd(i,j))/dy;
        else
            Cr_x(i,j)=0;
            Cr_y(i,j)=0;
        end

        a_x(i,j)=m*Cr_x(i,j)-fix(m*Cr_x(i,j));
        a_y(i,j)=m*Cr_y(i,j)-fix(m*Cr_y(i,j));

        a1_x(i,j)=a_x(i,j)^2*(3-2*a_x(i,j));
        a2_x(i,j)=1-a1_x(i,j);
        a3_x(i,j)=a_x(i,j)^2*(1-a_x(i,j))*dx;
        a4_x(i,j)=-a_x(i,j)*(1-a_x(i,j))^2*dx;

        a1_y(i,j)=a_y(i,j)^2*(3-2*a_y(i,j));
        a2_y(i,j)=1-a1_y(i,j);
        a3_y(i,j)=a_y(i,j)^2*(1-a_y(i,j))*dy;
        a4_y(i,j)=-a_y(i,j)*(1-a_y(i,j))^2*dy;

        b1_x(i,j)=6*a_x(i,j)*(a_x(i,j)-1)/dx;
        b2_x(i,j)=-b1_x(i,j);
        b3_x(i,j)=a_x(i,j)*(3*a_x(i,j)-2);
        b4_x(i,j)=(a_x(i,j)-1)*(3*a_x(i,j)-1);

        b1_y(i,j)=6*a_y(i,j)*(a_y(i,j)-1)/dy;
        b2_y(i,j)=-b1_y(i,j);
        b3_y(i,j)=a_y(i,j)*(3*a_y(i,j)-2);
        b4_y(i,j)=(a_y(i,j)-1)*(3*a_y(i,j)-1);
    end
end

figure(1);
contourf(x1,y1,c1')
xlabel("x");
ylabel("y");
title("Initial distribution");

figure(13);
mesh(c1');
xlabel("x");
ylabel("y");
title("Initial distribution");

% disp(max(max(seita)));
% disp(min(min(seita)));
ck=c1;
cm=c1;
cl=c1;
cn=c1;
cd=c1;
c2=c1;
cxk=cx1;
cxm=cx1;
cxd=cx1;
cyl=cy1;
cyn=cy1;
cyd=cy1;
cx2=cx1;
cy2=cy1;
cxy2=cxy1;
e1=c1;
ex=0;ey=3200;
er=zeros(Nx+1,Ny+1);

for i=1:Nx+1
    for j=1:Ny+1
        e1(i,j)=c1(i,j);
%         e1(i,j)=sx/sqrt(sx^2+2*D*T)*sy/(sqrt(sy^2+2*D*T))*exp(-(x1(i)-ex)^2/(2*(sx^2+2*D*T))-(y1(j)-ey)^2/(2*(sy^2+2*D*T)));
    end
end
% [X,Y]=meshgrid(1:Nx+1,1:Ny+1);%%%%meshgrid,注意矩阵行数对应y，列数对应X，比如3*4列，则y=1:3，x=1:4
figure(2)
contourf(x1,y1,c1')
% contourf(x1,y1,e1')
xlabel("x");
ylabel("y");
title("Exact solution");

figure(14)
plot(x1,c1(:,Ny/4*3+1),'b');
% figure(4)
% plot(c1(65,:),'b');

% figure(3);
% contourf(x1,y1,u');
% xlabel("x");
% ylabel("y");
% title("u");

% figure(4);
% contourf(x1,y1,v');
% xlabel("x");
% ylabel("y");
% title("v");

% figure(5);
% contourf(x1,y1,cx1');
% xlabel("x");
% ylabel("y");
% title("cx1");
% 
figure(11)
mesh(cx1');
xlabel("x");
ylabel("y");
title("cx1");
% 
% figure(6);
% contourf(x1,y1,cy1');
% xlabel("x");
% ylabel("y");
% title("cy1");
% 
% figure(7);
% contourf(x1,y1,cxy1');
% xlabel("x");
% ylabel("y");
% title("cxy1");
% 
% figure(13);
% mesh(cxy1');
% xlabel("x");
% ylabel("y");
% title("cxy1");

% figure(10);
% contourf(x1,y1,cxt');
% xlabel("x");
% ylabel("y");
% title("cxt");

% figure(12);
% mesh(cxt');


for n=1:M
    for i=3:Nx-1
        for j=3:Ny-1
%%%%%%%%%%  Find out the C, Cx, Cy value at point d
            if i>=Nx/2+1 && j>Ny/2+1 % for quadrant I & Y positive half shaft
                ck(i,j)=a1_y(i,j)*c1(i,j-1)+a2_y(i,j)*c1(i,j)+a3_y(i,j)*cy1(i,j-1)+a4_y(i,j)*cy1(i,j);
                cm(i,j)=a1_y(i,j)*c1(i+1,j-1)+a2_y(i,j)*c1(i+1,j)+a3_y(i,j)*cy1(i+1,j-1)+a4_y(i,j)*cy1(i+1,j);
                cxk(i,j)=a1_y(i,j)*cx1(i,j-1)+a2_y(i,j)*cx1(i,j)+a3_y(i,j)*cxy1(i,j-1)+a4_y(i,j)*cxy1(i,j);
                cxm(i,j)=a1_y(i,j)*cx1(i+1,j-1)+a2_y(i,j)*cx1(i+1,j)+a3_y(i,j)*cxy1(i+1,j-1)+a4_y(i,j)*cxy1(i+1,j);
                
                cl(i,j)=a1_x(i,j)*c1(i,j)+a2_x(i,j)*c1(i+1,j)+a3_x(i,j)*cx1(i,j)+a4_x(i,j)*cx1(i+1,j);
                cn(i,j)=a1_x(i,j)*c1(i,j-1)+a2_x(i,j)*c1(i+1,j-1)+a3_x(i,j)*cx1(i,j-1)+a4_x(i,j)*cx1(i+1,j-1);
                cyl(i,j)=a1_x(i,j)*cy1(i,j)+a2_x(i,j)*cy1(i+1,j)+a3_x(i,j)*cxy1(i,j)+a4_x(i,j)*cxy1(i+1,j);
                cyn(i,j)=a1_x(i,j)*cy1(i,j-1)+a2_x(i,j)*cy1(i+1,j-1)+a3_x(i,j)*cxy1(i,j-1)+a4_x(i,j)*cxy1(i+1,j-1);
            elseif i<Nx/2+1 && j>=Ny/2+1 % for quadrant II & X negative half shaft
                ck(i,j)=a1_y(i,j)*c1(i,j)+a2_y(i,j)*c1(i,j+1)+a3_y(i,j)*cy1(i,j)+a4_y(i,j)*cy1(i,j+1);
                cm(i,j)=a1_y(i,j)*c1(i+1,j)+a2_y(i,j)*c1(i+1,j+1)+a3_y(i,j)*cy1(i+1,j)+a4_y(i,j)*cy1(i+1,j+1);
                cxk(i,j)=a1_y(i,j)*cx1(i,j)+a2_y(i,j)*cx1(i,j+1)+a3_y(i,j)*cxy1(i,j)+a4_y(i,j)*cxy1(i,j+1);
                cxm(i,j)=a1_y(i,j)*cx1(i+1,j)+a2_y(i,j)*cx1(i+1,j+1)+a3_y(i,j)*cxy1(i+1,j)+a4_y(i,j)*cxy1(i+1,j+1);
                
                cl(i,j)=a1_x(i,j)*c1(i,j+1)+a2_x(i,j)*c1(i+1,j+1)+a3_x(i,j)*cx1(i,j+1)+a4_x(i,j)*cx1(i+1,j+1);
                cn(i,j)=a1_x(i,j)*c1(i,j)+a2_x(i,j)*c1(i+1,j)+a3_x(i,j)*cx1(i,j)+a4_x(i,j)*cx1(i+1,j);
                cyl(i,j)=a1_x(i,j)*cy1(i,j+1)+a2_x(i,j)*cy1(i+1,j+1)+a3_x(i,j)*cxy1(i,j+1)+a4_x(i,j)*cxy1(i+1,j+1);
                cyn(i,j)=a1_x(i,j)*cy1(i,j)+a2_x(i,j)*cy1(i+1,j)+a3_x(i,j)*cxy1(i,j)+a4_x(i,j)*cxy1(i+1,j);
            elseif i<=Nx/2+1 && j<Ny/2+1 % for quadrant III & Y negative half shaft
                ck(i,j)=a1_y(i,j)*c1(i-1,j)+a2_y(i,j)*c1(i-1,j+1)+a3_y(i,j)*cy1(i-1,j)+a4_y(i,j)*cy1(i-1,j+1);
                cm(i,j)=a1_y(i,j)*c1(i,j)+a2_y(i,j)*c1(i,j+1)+a3_y(i,j)*cy1(i,j)+a4_y(i,j)*cy1(i,j+1);
                cxk(i,j)=a1_y(i,j)*cx1(i-1,j)+a2_y(i,j)*cx1(i-1,j+1)+a3_y(i,j)*cxy1(i-1,j)+a4_y(i,j)*cxy1(i-1,j+1);
                cxm(i,j)=a1_y(i,j)*cx1(i,j)+a2_y(i,j)*cx1(i,j+1)+a3_y(i,j)*cxy1(i,j)+a4_y(i,j)*cxy1(i,j+1);
                
                cl(i,j)=a1_x(i,j)*c1(i-1,j+1)+a2_x(i,j)*c1(i,j+1)+a3_x(i,j)*cx1(i-1,j+1)+a4_x(i,j)*cx1(i,j+1);
                cn(i,j)=a1_x(i,j)*c1(i-1,j)+a2_x(i,j)*c1(i,j)+a3_x(i,j)*cx1(i-1,j)+a4_x(i,j)*cx1(i,j);
                cyl(i,j)=a1_x(i,j)*cy1(i-1,j+1)+a2_x(i,j)*cy1(i,j+1)+a3_x(i,j)*cxy1(i-1,j+1)+a4_x(i,j)*cxy1(i,j+1);
                cyn(i,j)=a1_x(i,j)*cy1(i-1,j)+a2_x(i,j)*cy1(i,j)+a3_x(i,j)*cxy1(i-1,j)+a4_x(i,j)*cxy1(i,j);
                
            elseif i>Nx/2+1 && j<=Ny/2+1 % for quadrant IV & X positive half shaft
                ck(i,j)=a1_y(i,j)*c1(i-1,j-1)+a2_y(i,j)*c1(i-1,j)+a3_y(i,j)*cy1(i-1,j-1)+a4_y(i,j)*cy1(i-1,j);
                cm(i,j)=a1_y(i,j)*c1(i,j-1)+a2_y(i,j)*c1(i,j)+a3_y(i,j)*cy1(i,j-1)+a4_y(i,j)*cy1(i,j);
                cxk(i,j)=a1_y(i,j)*cx1(i-1,j-1)+a2_y(i,j)*cx1(i-1,j)+a3_y(i,j)*cxy1(i-1,j-1)+a4_y(i,j)*cxy1(i-1,j);
                cxm(i,j)=a1_y(i,j)*cx1(i,j-1)+a2_y(i,j)*cx1(i,j)+a3_y(i,j)*cxy1(i,j-1)+a4_y(i,j)*cxy1(i,j);
                
                cl(i,j)=a1_x(i,j)*c1(i-1,j)+a2_x(i,j)*c1(i,j)+a3_x(i,j)*cx1(i-1,j)+a4_x(i,j)*cx1(i,j);
                cn(i,j)=a1_x(i,j)*c1(i-1,j-1)+a2_x(i,j)*c1(i,j-1)+a3_x(i,j)*cx1(i-1,j-1)+a4_x(i,j)*cx1(i,j-1);
                cyl(i,j)=a1_x(i,j)*cy1(i-1,j)+a2_x(i,j)*cy1(i,j)+a3_x(i,j)*cxy1(i-1,j)+a4_x(i,j)*cxy1(i,j);
                cyn(i,j)=a1_x(i,j)*cy1(i-1,j-1)+a2_x(i,j)*cy1(i,j-1)+a3_x(i,j)*cxy1(i-1,j-1)+a4_x(i,j)*cxy1(i,j-1);
            end
            
%             cd(i,j)=a1_x(i,j)*ck(i,j)+a2_x(i,j)*cm(i,j)+a3_x(i,j)*cxk(i,j)+a4_x(i,j)*cxm(i,j);
%             cd(i,j)=a1_y(i,j)*cn(i,j)+a2_y(i,j)*cl(i,j)+a3_y(i,j)*cyn(i,j)+a4_y(i,j)*cyl(i,j);
            cd(i,j)=1/2*(a1_x(i,j)*ck(i,j)+a2_x(i,j)*cm(i,j)+a3_x(i,j)*cxk(i,j)+a4_x(i,j)*cxm(i,j))+1/2*(a1_y(i,j)*cn(i,j)+a2_y(i,j)*cl(i,j)+a3_y(i,j)*cyn(i,j)+a4_y(i,j)*cyl(i,j));
            cxd(i,j)=b1_x(i,j)*ck(i,j)+b2_x(i,j)*cm(i,j)+b3_x(i,j)*cxk(i,j)+b4_x(i,j)*cxm(i,j);
            cyd(i,j)=b1_y(i,j)*cn(i,j)+b2_y(i,j)*cl(i,j)+b3_y(i,j)*cyn(i,j)+b4_y(i,j)*cyl(i,j);
        end
    end
%%%%%%%%%%  Pass the C value of point d to node(i,j)
%     for i=3:Nx-1
%         for j=3:Ny-1
%             c1(i,j)=cd(i,j);
%         end
%     end
    c1=cd;
    
%     diffu_x=c1(i+1,j)-2*c1(i,j)+c1(i-1,j);
%     diffu_y=c1(i,j+1)-2*c1(i,j)+c1(i,j-1);
%     for i=2:Nx
%         for j=2:Ny
%             c1(i,j)=c1(i,j)+(dt/(dx^2))*D*200*diffu_x+(dt/(dy^2))*D*diffu_y*200;
%         end
%     end

%     if mod(M,300)==0
%         for i=3:Nx-1
%             for j=3:Ny-1
%                 if c1(i,j)>0.1
%                     show(i,j)=c1(i,j);
%                 end
%             end
%         end
%     end
        
%%%%%%%%%%  Calculate the Cx, Cy, Cxy value of node(i,j) with the value at
%%%%%%%%%%  point d
    for i=3:Nx-1
        for j=3:Ny-1
           cx2(i,j)=cxd(i,j)-w*dt*1/2*(cy1(i,j)+cyd(i,j));
           cy2(i,j)=cyd(i,j)+w*dt*1/2*(cx1(i,j)+cxd(i,j));
%            cx2(i,j)=cxd(i,j)-w*dt*cyd(i,j);
%            cy2(i,j)=cyd(i,j)+w*dt*cxd(i,j);
        end
    end
    cx1=cx2;cy1=cy2;
    
    for i=3:Nx-1
        for j=3:Ny-1
           cxy2(i,j)=1/2*(cx1(i,j+1)-cx1(i,j-1))/(y1(j+1)-y1(j-1)) + 1/2*(cy1(i+1,j)-cy1(i-1,j))/(x1(i+1)-x1(i-1));
        end
    end
    cxy1=cxy2;
    
%     if mod(n,10)==0
%         figure(10);
%         hold on;
%         xlabel("x");
%         ylabel("y");
%         contourf(x1,y1,c1');
%         title(n+"/"+M);
%     end
    figure(10);
    xlabel("x");
    ylabel("y");
    hold on
    if mod(n,20)==0
        contourf(x1,y1,c1');
        title(n+"/"+M);
    end
    hold off
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,128);
    if n==1
         imwrite(imind,cm,'square1.gif','gif', 'Loopcount',inf,'DelayTime',0.001);
    end
    if mod(n,5)==0
         imwrite(imind,cm,'square1.gif','gif','WriteMode','append','DelayTime',0.001);
    end
end

toc;
% figure(1)
% hold on;
% plot(c1(:,64),'r');

figure(14)
hold on
% plot(x1,e1(:,Ny/4*3+1),'g');
plot(x1,c1(:,Ny/4*3+1),'r');
xlabel("x");
ylabel("C");
% legend({"Initial distribution","Exact solution",Nx+"x"+Ny},"Location","northwest");
legend({"Exact solution",Nx+"x"+Ny},"Location","northwest");
title(Nx+"x"+Ny);
maxvalue=max(max(c1));
text(-7500,0.75,"max:"+maxvalue);
% text(-7500,0.75,"max:0.9990");
hold off

figure(8);
contourf(x1,y1,c1');
xlabel("x");
ylabel("y");
title(Nx+"x"+Ny);

figure(9);
xlabel("x");
ylabel("y");
title("C")
mesh(x1,y1,c1');

% figure(15);
% contourf(x1,y1,show');
% xlabel("x");
% ylabel("y");
% title(Nx+"x"+Ny);


error = 0;
for i=1:Nx+1
    for j=1:Ny+1
        error = error + abs(c1(i,j)-e1(i,j))^2;
    end
end
error = error/((Nx+1)*(Ny+1));
disp(error);
disp(log(error));
disp(sqrt(dx^2+dy^2));
disp(log(sqrt(dx^2+dy^2)));

% x=[-3.8123 -4.5055 -5.1986 -5.8918];
% y=[-5.4288 -6.6794 -7.8829 -9.1738];
% figure(15);
% plot(x,y,'-o');
% xlabel("log(sqrt(dx^2+dy^2))");
% ylabel("log(Error)");
% %最简单的线性回归matlab代码
% p=polyfit(x,y,1)%1代表希望拟合出的曲线是线性的
% 
% p =
% 
%     1.7945    1.4156
