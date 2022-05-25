% Supercell edge modes
% clear all
% Number of Resonators
k = 4;
ll = 1;
N = 6*k*ll;
% N_modulated = ceil(k/2)*6;
N_modulated = N;

phase_shift = [repmat([2/3,4/3,0,2/3,0,4/3]*pi,1,k*ll)];    

Omega = 2;
eps = 0;
dd = 0;
[V_mode_0,f_exp_temp,f_exp_im,ind,c,M_0] = produce_edge_supercell_location(k,ll,dd,Omega,phase_shift);

figure
for i = 1:N
    subplot(3,round(N/3)+1,i)
    plot(1:N,real(V_mode_0(1:N,i)));
end
close

edge_ind = [k*2 k*2+1];

figure
plot(1:N,real(V_mode_0(1:N,edge_ind(1))),'.k', 'MarkerSize',20);
xlabel('The serial number of the resonators')
ylabel('Re(\Psi)')
figure
plot(1:N,real(V_mode_0(1:N,edge_ind(2))),'.k', 'MarkerSize',20);
xlabel('The serial number of the resonators')
ylabel('Re(\Psi)')
% 
function [V,f_exp,f_exp_im,ind,c,M_0] = produce_edge_supercell_location(k,ll,dd,Omega,phase_shift)
N = 6*k*ll;

% radii
    R = 0.1*ones(1,N);

    %%% Material parameters
    rho0 = 1e3;             % background material
    kappa0 = 1e3;           % background material
    v = sqrt(kappa0/rho0);  % speed of sound in water

    rho_b = 1;            % density of resonators  
    kappa_b = 1;          % bulk modulus of resonators
    v_b = sqrt(kappa_b/rho_b);  % speed of sound resonators

    
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
c_unit = [c1+l(2*pi/3)+dd;c1+l(4*pi/3);c1+l(0);c2+l(pi/3);c2+l(3*pi/3);c2+l(5*pi/3)];
L_unit = sqrt(3);
     c = c_unit;
    for ii = 1:k-1
        c = [c ; c_unit + ii*[L_unit 0 ; L_unit 0;L_unit 0;L_unit 0 ; L_unit 0;L_unit 0]];
    end


        
    c = sortrows(c,1);
    cx = c(:,1);
    cy = c(:,2);
    
    % Plot the geometry
%     [cy,indc]=sort(c(:,2));
%     c(:,1)=c(indc,1);
%     c(:,2)=cy;
%     cx = c(:,1);
%     cy = c(:,2);
    figure, hold on
    t = linspace(0,2*pi);
    for n = 1:N
        plot(cx(n)+R(n)*cos(t), cy(n)+R(n)*sin(t),'k')
        text (cx(n), cy(n), num2str(n))
    end
    
    daspect([1 1 1])
    hold off
%     close
    
    eps_modulation = zeros(1,N);
    rhot = @(t) rho_b./(1 + eps_modulation.*cos(Omega*t + phase_shift)); 

    sqrtkappat = @(t) sqrt(kappa_b)*ones(N,1); % ./[1 + eps*sin(t + pi*(1:N))];
    dkappa = @(t) 0*R;
    d2kappa = @(t) 0*R;
    w3 = @(t) 0*(1:N);

    T = 2*pi/Omega;
    steps = 1000;

    % High contrast parameters \delta
    delta=rho_b/rho0;


 
    %% Compute the resonances using the capacitance matrix, which is approximated using the multipole expansion method
    C = capacitance_2D(c,R,rho0,rho_b,kappa0,kappa_b,delta);
%     vol = pi*R.^2;
    vol = 4/3*pi*R.^3;
    vol = vol';
    CoeffMat = @(t) makeM(t,delta,kappa0,rho0,vol,C,rhot, sqrtkappat, w3);
    M_0 = delta*kappa0/rho0*inv(diag(vol))*C;
    %% Solve for Psi
    [Tspan, X_fundamental] = HillSolver(CoeffMat,T,steps); % X_fundamental = [Psi; dPsi/dt]
    [TOUT, X_bloch,V_mode,V,D] = Hill2BlochModes(CoeffMat,T,steps);
    
    
    %% Solve Floquet exp
    X_T = squeeze(X_fundamental(end,:,:));
    [V,d] = eig(X_T,'vector');
    [f_exp,ind] = sort(imag(log(d)/T),'descend');
    f_exp_im = real(log(d(ind))/T);
    V = V(:,ind);
end