%================================
%
%  plot data from HYPRE test
%
%================================


function HyprePlot

clear all
close all
clc

% read in psi, rhs, u
testCase = 1;
if (testCase == 0)
    u   = importdata('exJanaSquare_u.dat');
    psi = importdata('exJanaSquare_psi.dat');
    rhs = importdata('exJanaSquare_rhs.dat');
else
    u   = importdata('exJanaSphere_u.dat');
    psi = importdata('exJanaSphere_psi.dat');
    rhs = importdata('exJanaSphere_rhs.dat');
end

[nx,ny] = size(u);


% get analytical solution
ax = 0; ay = 0;
bx = 1; by = 1;
CCgrid = 0;
[X,Y,h] = initGrid(nx,ax,bx,ay,by,CCgrid);
[uexact, fAna] = getExactSolution(X,Y,testCase);

nx
ny
size(X)


%========================================
% plot output

% figure set up:
fs  = 15;
fid = 1;

plot_cut = 1;
plot_psi = 1;
plot_rhs = 1;

%------------ numerical solution --------------------
figure('color','w')
figure(fid)
set(gca,'Fontsize',fs);
pcolor(u.*psi);
shading interp;
colormap default
title({'u num';['Nx = ',num2str(nx)]})
xlabel('x')
ylabel('y')
fid=fid+1;
colorbar

%------------ Uex --------------------
figure('color','w')
figure(fid)
set(gca,'Fontsize',fs);
pcolor(uexact.*psi);
shading interp;
colormap default
title({'Uexact \Omega_1';['Nx = ',num2str(nx)]})
xlabel('x')
ylabel('y')
fid=fid+1;
colorbar


if (plot_cut)
    mid = round(nx*0.5);
   figure('color','w')
   figure(fid)
   hold on
   set(gca,'Fontsize',fs); 
   plot(X(:,1),uexact(:,mid).*psi(:,mid),'-r*')
   plot(X(:,1),u(:,mid).*psi(:,mid),'-bo')
   legend('exact','num.')
   xlabel('x');
   ylabel('u(x)')
   title('Cut through middle')
   box on; grid on;
   fid=fid+1;
   
      figure('color','w')
   figure(fid)
   set(gca,'Fontsize',fs); 
   plot(X(:,1),psi(:,mid),'-r*')
   xlabel('x');
   ylabel('\psi(x)')
   title('PSI cut')
   box on; grid on;
   fid=fid+1;
   
   
end;

    

if (plot_psi)
    %------------ psi --------------------
    figure('color','w')
    figure(fid)
    set(gca,'Fontsize',fs);
    pcolor(psi);
    shading interp;
    colormap default
    title({'Phase Field in \Omega_1';['Nx = ',num2str(nx)]})
    xlabel('x')
    ylabel('y')
    fid=fid+1;
    colorbar
end;

if (plot_rhs)
    
    %------------ f anal --------------------
    h2 = h*h;
    figure('color','w')
    figure(fid)
    set(gca,'Fontsize',fs);
    pcolor(fAna.*psi.*h2);
    shading interp;
    colormap default
    title({'Fana \Omega_1';['Nx = ',num2str(nx)]})
    xlabel('x')
    ylabel('y')
    fid=fid+1;
    colorbar
    
    %------------ rhs num --------------------
    figure('color','w')
    figure(fid)
    set(gca,'Fontsize',fs);
    pcolor(rhs);
    shading interp;
    colormap default
    title({'RHS in \Omega_1';['Nx = ',num2str(nx)]})
    xlabel('x')
    ylabel('y')
    fid=fid+1;
    colorbar
end;