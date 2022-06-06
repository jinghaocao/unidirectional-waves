% Mirrored Supercell edge modes, used to produce Figure 7 and 8.

% Number of Resonators
k = 5;
ll = 1;
N = 12*k*ll;


% Perturbed case, change eps = 0 to obatian the unperturbed case
eps = 0.4;
Omega = 2;
[V,f_exp,f_exp_im,ind,c,M_0] = produce_edge_supercell(k,ll,eps,Omega);

% Plot all the bloch modes
figure 
hold on 
for j = 1:N
    subplot(5,ceil(N/5),j)
    plot(1:N,real(V(1:N,j)));
    xlabel(strcat('The serial number of the resonators',' ,eps=',num2str(eps)))
    ylabel('Re(\Psi)')
end

% Define the edge index and plot all the bloch edge modes
edge_ind = [1 k*4 k*4+1];
kk = 1;
figure 
hold on 
for j = edge_ind
    subplot(3,1,kk)
    plot(1:N,real(V(1:N,j)));
    xlabel(strcat('The serial number of the resonators',' ,eps=',num2str(eps)))
    ylabel('Re(\Psi)')
    kk = kk+1;
end


function [V,f_exp,f_exp_im,ind,c,M_0] = produce_edge_supercell(k,ll,eps,Omega)
    N = 6*k*ll*2;
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
    d = [c2(1)-c1(1) 0];
    eps = 3*R(1); %eps = 3*R;
    l = @(theta) eps*[cos(theta+pi/6),sin(theta+pi/6)];
    c_unit = [c1+l(2*pi/3);c1+l(4*pi/3);c1+l(0);c2+l(pi/3);c2+l(3*pi/3);c2+l(5*pi/3)];
    c_unit_2 = [c1+d+l(2*pi/3);c1+d+l(4*pi/3);c1+d+l(0);c2-d+l(pi/3);c2-d+l(3*pi/3);c2-d+l(5*pi/3)];

    L_unit = sqrt(3);
    c = c_unit;
    for ii = 1:k-1
        c = [c ; c_unit + ii*[L_unit 0 ; L_unit 0;L_unit 0;L_unit 0 ; L_unit 0;L_unit 0]];
    end

    for ii = k:2*k-1
    c = [c ; c_unit_2 + ii*[L_unit 0 ; L_unit 0;L_unit 0;L_unit 0 ; L_unit 0;L_unit 0]];
    end

    for ii = 1:ll-1
        c = [c ; c + ii*repmat([0 1],k*6,1)];
    end

    phase_shift_left = [repmat([2/3,4/3,0,2/3,0,4/3]*pi,1,k*ll)];    
    phase_shift_right = [repmat([4/3,0,2/3,0,4/3,2/3]*pi,1,k*ll)];
    phase_shift = [phase_shift_left phase_shift_right];
    c = sortrows(c,1);
    cx = c(:,1);
    cy = c(:,2);
    
    % Plot the geometry
    figure, hold on
    t = linspace(0,2*pi);
    for n = 1:N
        plot(cx(n)+R(n)*cos(t), cy(n)+R(n)*sin(t),'k')
        text (cx(n), cy(n), num2str(n))
    end
    
    daspect([1 1 1])
    hold off
    close

    eps_modulation = eps*ones(1,N);
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