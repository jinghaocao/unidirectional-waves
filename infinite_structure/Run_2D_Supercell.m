


%%% Computes the quasifrequencies for resonator-modulated systems with
%%% hexagonal lattice.

%close all;

%%% lattice vectors ( a1 = [L1x;0], a2 = [L2x;L2y], D: distance between nearest neighbor bubbles ) 
D = 1; 
L1x = D*sqrt(3);
L2x = D*sqrt(3)/2;
L2y = D*3/2;
L1 = [L1x,0];
L2 = [L2x,L2y];

%%% Set reciprocal vectors
b1 = 4*pi/(3*D)*[sqrt(3)/2, -1/2];
b2 = 4*pi/(3*D)*[0, 1];

%%% Set symmetry points in reciprocal space
pointG = [0,0];
pointM = 1/2*b1;
pointK = 2/3*b1 + 1/3*b2;

%%% Radius of bubble
R=D*0.1; %R = D*0.1;
vol = pi*R^2;

%%% Material parameters for bubble
rho_0 = 9000;
rho_b = 1;
kappa_0 = 9000;
kappa_b = 1;

%%% High contrast parameter \delta
delta = rho_b/rho_0;
weight = delta*kappa_b/rho_b/vol;

%%% Set positions of two bubbles in unit cell



N = 6; % Number of resonators in unit cell
    % supercell
    D = 1; 
L1x = D*sqrt(3);
L2x = D*sqrt(3)/2;
L2y = D*3/2;
L1 = [L1x,0];
L2 = [L2x,L2y];
c1 = 1/3*(L1+L2);
c2 = 2/3*(L1+L2);
eps = 3*R(1); %eps = 3*R;
l = @(theta) eps*[cos(theta+pi/6),sin(theta+pi/6)];
c = [c1+l(2*pi/3)+dd;c1+l(4*pi/3);c1+l(0);c2+l(pi/3);c2+l(3*pi/3);c2+l(5*pi/3)];

%%% Plot the geometry
figure, hold on
t = linspace(0,2*pi);
for i = -2:2
    for j = -2:2
        for n = 1:N
            plot(c(n,1)+R*cos(t)+i*L1(1)+j*L2(1), c(n,2)+R*sin(t)+i*L1(2)+j*L2(2),'k')  
            text(c(n,1),c(n,2),num2str(n))
        end
    end
end
daspect([1 1 1])
hold off

%%% Define the time modulation

%%% Cones using R = 0.02*D, eps = 13*R
% eps = 0 ; 
% epsk = 0.51;
% Omega = 1.5;

% eps = 0;
% epsk = 0.237;
% Omega = 0.5;
% 
% eps = 0.38;
% epsk = 0;
% Omega = 0.5;

%%% Cone using R=0.1*D, eps = 3*R 
% eps = 0.3;
% epsk = 0;
% Omega = 0.15;


eps = 0; % eps = 0.2;
epsk = 0;
Omega = 0.15;

modr = @(t,phi) 1./(1+eps*cos(Omega*t+phi));
rhot = @(t) [modr(t,0); modr(t,2*pi/3); modr(t,4*pi/3); modr(t,0); modr(t,2*pi/3); modr(t,4*pi/3)];
modk = @(t,phi) sqrt(1+epsk*cos(Omega*t+phi));
w3k = @(t,phi) epsk*Omega^2*(-4*cos(Omega*t+phi) +epsk*(-5+cos(2*(Omega*t+phi))) ) / (8*(1 + epsk*cos(Omega*t+phi))^2);
%w3k = @(t,phi) -3/4*(epsk*cos(Omega*t+phi)/(1+epsk*sin(Omega*t+phi))).^2 + epsk*sin(Omega*t+phi)/(1+epsk*sin(Omega*t+phi));
sqrtkappat = @(t) [modk(t,0); modk(t,2*pi/3); modk(t,4*pi/3); modk(t,0); modk(t,2*pi/3); modk(t,4*pi/3)];
w3 = @(t)  [w3k(t,0); w3k(t,2*pi/3); w3k(t,4*pi/3); w3k(t,0); w3k(t,2*pi/3); w3k(t,4*pi/3)];
T = 2*pi/Omega;

%%% Set the discretization or truncation parameters
N_multipole=3;
N_lattice=30;

N_KM = 80;
N_MG = floor(sqrt(3)*N_KM);
N_GK = 2*N_KM;

%%% Precompute values of zeta function (used in computation of lattice sums, just for acceleration)
d_zeta=makezetadata;
BandDataKM = zeros(2*N,N_KM);
BandDataMG = zeros(2*N,N_MG);
BandDataGK = zeros(2*N,N_GK);

k0 = 0.001;
JHdata = makeJHdata0original(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);
    
parfor alp_ind = 1:N_MG
    %%% Set quasi-periodic parameter \alpha
    alp =  pointM + (alp_ind-1)/N_MG * (pointG-pointM);

    %%% Compute the capacitance matrix C
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    BandDataMG(:, alp_ind) = hill_exp(T,M,N);
end

parfor alp_ind = 1:N_GK    
    %%% Set quasi-periodic parameter \alpha
    alp =  pointG + (alp_ind-1)/N_GK * (pointK-pointG);

    %%% Compute the operator A
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    BandDataGK(:, alp_ind) = hill_exp(T,M,N);
end

parfor alp_ind = 1:N_KM
    %%% Set quasi-periodic parameter \alpha
    alp =  pointK + (alp_ind-1)/N_KM * (pointM-pointK);

    %%% Compute the operator A
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    BandDataKM(:, alp_ind) = hill_exp(T,M,N);
end    


save BandDataMG
save BandDataGK
save BandDataKM

lengthMG = norm(pointG-pointM);
lengthGK = norm(pointK-pointG);
lengthKM = norm(pointM-pointK);
lengthMGK = lengthMG + lengthGK + lengthKM;
alps = [(0:N_MG-1)*lengthMG/N_MG, (0:N_GK-1)*lengthGK/N_GK + lengthMG, (0:N_KM-1)*lengthKM/N_KM + lengthMG + lengthGK];

%%

% figure
% hold on
% %w_end=ceil(10*max(max(abs([BandDataMG, BandDataGK, BandDataKM]))))/10;
% %axis([0 lengthMGK 0 w_end])
% for j=1:N_MG
%     for i = 1:2*N
%         X=(j-1)*lengthMG/N_MG;
%         Y=BandDataMG(i,j);
%         if real(Y)>0
%         plot(X,real(Y),'.','color',[0.2 0.3 0.7],'MarkerSize',5)
%         plot(X,imag(Y),'r.','MarkerSize',5)
%         end
%     end
% end
% 
% for j=1:N_GK
%     for i=1:2*N
%         X=(j-1)*lengthGK/N_GK + lengthMG;
%         Y=BandDataGK(i,j);
%         if real(Y)>0
%         plot(X,real(Y),'.','color',[0.2 0.3 0.7],'MarkerSize',5)
%         plot(X,imag(Y),'r.','MarkerSize',5)
%         end
%     end
% end
% 
% for j=1:N_KM
%     for i=1:2*N
%         X=(j-1)*lengthKM/N_KM + lengthMG + lengthGK;
%         Y=BandDataKM(i,j);
%         if real(Y)>0
%         plot(X,real(Y),'.','color',[0.2 0.3 0.7],'MarkerSize',5)
%         plot(X,imag(Y),'r.','MarkerSize',5)
%         end
%     end
% end
% 
% hold on
% x = [ 0*pi, lengthMG ,  lengthMG+lengthGK,      lengthMGK];
% y= ['    M ';'\Gamma';'   K  ';'  M   '];
% 
% set(gca,'XTick',x); % Change x-axis ticks
% set(gca,'XTickLabel',y); 
% ylabel('\omega')
% %title('Sub-wavelength bands for the honeycomb bubble structure','FontWeight','normal')
% 
% set(gca,'fontsize', 15);
% set(gca,'xgrid','on')

%% 
a_plot = [(0:N_MG-1)*lengthMG/N_MG, (0:N_GK-1)*lengthGK/N_GK + lengthMG, (0:N_KM-1)*lengthKM/N_KM + lengthMG + lengthGK];
w_plot = [BandDataMG, BandDataGK, BandDataKM];
N_a = length(a_plot);


% Plots the solution and saves good images
fsz = 16;
lfsz = 14;
alw = 0.75;
lw = 1.2;

figure
hold on
w_real = sort(real(w_plot),1);
w_imag = sort(imag(w_plot),1);

for I_a=1:N_a
    for I_w=N+1:2*N
        if w_real(I_w,I_a) > Omega/2*(1-1e-6)
            w_real(I_w,I_a) = w_real(I_w,I_a) - Omega;
        end
    end
end
%w_real(2:end,292) = w_real(1:end-1,292); w_real(1,292) = - w_real(end,292); %only at eps = 0.2
%w_real(2:end,306) = w_real(1:end-1,306); w_real(1,306) = - w_real(end,306);

for I_w = 1:2*N
plot(a_plot,w_real(I_w,:),'-b','LineWidth',lw)
%plot(a_plot,w_imag(I_w,:),'--r','LineWidth',lw)
end
x = [ 0*pi, lengthMG ,  lengthMG+lengthGK,      lengthMGK];
y= ['     M  ';'$\Gamma$';'    K   ';'   M    '];

axis([0,lengthMGK,0,Omega/2])
set(gca,'XTick',x); % Change x-axis ticks
set(gca,'XTickLabel',y); 

set(gca,'fontsize', 15);
set(gca,'xgrid','on')
%legend('Real part','Imaginary part','position',[0.73,0.63,0.2,0.1],'interpreter','latex');
xlabel('Quasi-periodicity $\alpha$','interpreter','latex');
ylabel('Frequency $\omega$','interpreter','latex');
ax = gca;

set(gcf, 'Position', [0, 0, 600, 500])
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(gca,'TickLabelInterpreter','latex')
%l.FontSize = lfsz;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height*0.99];
%print('HLresmod_b','-depsc');

%axis([lengthMG*0.5,lengthMG*1.5, 0, 0.05])
%print('HLresmod_c','-depsc');

% 4.40699 4.60348
%axis([0 lengthMGK 0.5 0.6])
%%% Assemble matrix M
function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
R = diag(rhot(t));
Rinv = diag(1./rhot(t));
K = diag(sqrtkappat(t));
W3 = diag(w3(t));
out = delta*kappab/(vol*rhob)*K*R*C*K*Rinv + W3;
end