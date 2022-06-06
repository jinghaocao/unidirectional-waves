% Computes the band functions of the 1D infinite structure 
% of supercell, this produces the figure 6 of the paper

% Define the modulation parameters
Omega = 0.2; 

L_unit = sqrt(3);

%%% Set reciprocal vectors
b1 = pi/L_unit;

%%% Set symmetry points in reciprocal space
pointG = 0;
pointM = b1;
N = 6;
epsr = 0.3*ones(1,N);

%%% Radius of bubble
R=0.1*ones(N,1);
vol = pi*R.^2;
%%% Material parameters for bubble
rho_0 = 9000;
rho_b = 1;
kappa_0 = 9000;
kappa_b = 1;

rho0 = 1e3;             % background material
kappa0 = 1e3;           % background material
    
%%% High contrast parameter \delta
delta = rho_b/rho_0;
% weight = delta*kappa_b/rho_b/vol;

%%% Set positions of two bubbles in unit cell
 % supercell
N=6;
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
c = [c1+l(2*pi/3);c1+l(4*pi/3);c1+l(0);c2+l(pi/3);c2+l(3*pi/3);c2+l(5*pi/3)];%% Plot the geometry
c = sortrows(c,1);

figure, hold on
t = linspace(0,2*pi);
for i = 0:2
    for n = 1:N
        plot(c(n,1)+R(n)*cos(t)+i*L_unit, c(n,2)+R(n)*sin(t),'k')  
        text(c(n,1),c(n,2),num2str(n));
    end
end

daspect([1 1 1])
hold off
close

%%% Define the time modulation
phii = [2/3,4/3,0,2/3,0,4/3]*pi;
rhot = @(t) 1./(1+epsr.*cos(Omega*t+phii));

epsk= 0;
sqrtkappat = @(t) [1./sqrt(1+cos(Omega*t+[0:1/(N-1):1]*pi).*epsk)];
w3 = @(t) [Omega^2*(4*cos(Omega*t+[0:1/(N-1):1]*0.5*pi).*epsk+(3+cos(2*(Omega*t+[0:1/(N-1):1]*0.5*pi))).*epsk)./(8*(1+cos(Omega*t+[0:1/(N-1):1].*0.5*pi).*epsk).^2)];
T = 2*pi/Omega;

%%% Set the discretization or truncation parameters
N_multipole=3;
N_a = 100;
% alps = linspace(-pointM,pointM,N_a);
alps = [linspace(-pointM,pointM,N_a)];

%%% Precompute values of zeta function (used in computation of lattice sums, just for acceleration)
d_zeta=makezetadata;
BandData = zeros(2*N,N_a);
cnd = zeros(1,N_a);

k0 = 0.001;
[Jdata,Hdata] = makeJHdata0(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);
for alp_ind = 1:N_a
    %%% Set quasi-periodic parameter \alpha
    alp =  alps(alp_ind);
    %%% Compute the operator A
    C = makeC_1D(k0,R,alp,L_unit,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);

    [BandData(:, alp_ind), cnd(alp_ind)] = hill_exp(T,M,N);
end  

%% Plots the solution and saves good images
fsz = 16;
lfsz = 14;
alw = 0.75;
lw = 1.2;

[w_real_temp, I_sort] = sort(real(BandData)); % sort(real(bands),1);
for I_a=1:N_a
    for I_w=1:2*N
        w_imag_temp(I_w,I_a) = imag(BandData(I_sort(I_w,I_a),I_a)); % sort(imag(bands),1);
    end
end

Ias = find(max(abs(w_imag_temp))>1e-4);
[w_i_sort, I_sort] = sort(w_imag_temp(N+1:2*N,Ias));

w_imag = sort(w_imag_temp);
w_real = w_real_temp;

for I_a=1:N_a
    for I_w=1:N*2
        if w_real(I_w,I_a) > Omega/2*(1-1e-6)
            w_real(I_w,I_a) = w_real(I_w,I_a) - Omega;
        end
    end
end
figure
hold on
for I_w = N+1:2*N
plot(alps,w_real(I_w,:),'.b','LineWidth',lw)
% plot(alps,w_imag(I_w,:),'-.r','LineWidth',lw)
end

% PLot band gap
% 
for j = N+1:2*N-2
plot(alps(2:N_a/2),max(w_real(j,2:N_a/2))*ones(1,(N_a/2-1)),'r')
plot(alps(2:N_a/2),min(w_real(j+1,2:N_a/2))*ones(1,(N_a/2-1)),'r')
plot(alps(N_a/2+1:N_a-1),max(w_real(j,N_a/2+1:N_a-1))*ones(1,(N_a/2-1)),'r')
plot(alps(N_a/2+1:N_a-1),min(w_real(j+1,N_a/2+1:N_a-1))*ones(1,(N_a/2-1)),'r')
end
% ylim([min([w_real(7,2),w_real(8,2),w_real(7,N_a-1),w_real(8,N_a-1)]),max([w_real(7,2),w_real(8,2),w_real(7,N_a-1),w_real(8,N_a-1)])]);

x = [ -pointM, 0,  pointM];
y= ['$-\pi/L$';'$\Gamma$';' $\pi/L$'];

xlim([-pointM,pointM])

set(gca,'XTick',x); % Change x-axis ticks
set(gca,'XTickLabel',y); 

set(gca,'xgrid','on')

%legend('Real part','Imaginary part','position',[0.74,0.62,0.2,0.1],'interpreter','latex','FontSize',lfsz);
xlabel('Quasiperiodicity $\alpha$','interpreter','latex');
ylabel('Frequency $\omega$','interpreter','latex');
ax = gca;

set(gcf, 'Position', [0, 0, 600, 500])
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(gca,'TickLabelInterpreter','latex')
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left*1.15 bottom ax_width*0.99 ax_height];
%print('SLresmod_b','-depsc');
print('oneDunfolded','-depsc');

%%% Assemble matrix M
function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
Rho = diag(rhot(t));
Rinv = diag(1./rhot(t));
K = diag(sqrtkappat(t));
W3 = diag(w3(t));
out = delta*kappab/rhob*inv(diag(vol))*K*Rho*C*K*Rinv + W3;
end