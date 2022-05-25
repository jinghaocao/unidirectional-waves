k=10;
N = 4*k+1;
figure
hold on
for eps = [0:0.1:0.3]
    eps_modulation = [ones(1,N)*eps];
    [f_exp,V] = produce_fexp_edge(N,eps_modulation);
    plot(1:N,real(V(1:N,N-2)))  
    xlim([1,N])
end
legend('\epsilon = 0','\epsilon = 0.1','\epsilon = 0.2','\epsilon = 0.3')

% eps_modulation = [ones(1,ceil(N/2))*eps, zeros(1,floor(N/2))];
% eps_modulation = [zeros(1,N)];
% 
% [f_exp,V] = produce_fexp_edge(N,eps_modulation);
% figure
% for j = 1:N
%     subplot(5,ceil(N/5),j)
%     plot(1:N,real(V(1:N,j)))
%     xlim([1,N])
% end



function [f_exp,V] = produce_fexp_edge(N,eps_modulation)
    % Number of Resonators

    % radii
    R = 0.1*ones(1,N);
    Omega = 0.2;

    %%% Material parameters
    rho0 = 1e3;             % background material
    kappa0 = 1e3;           % background material
    v = sqrt(kappa0/rho0);  % speed of sound in water

    rho_b = 1;            % density of resonators  
    kappa_b = 1;          % bulk modulus of resonators
    v_b = sqrt(kappa_b/rho_b);  % speed of sound resonators

    d1 = 1;
    d2 = 0.6;
    k = (N-1)/4;
    c_1 = [0 0.4];
    cx = [];
    for j = 0:k-1
      cx= [cx c_1+j*d1];  
    end
    mid = cx(end)+d2;
    cx = [cx mid];
    for j = 0:k-1
        cx = [cx c_1+mid+d2+j*d1];
    end
    phase_shift = [];
    k_left = floor(N/(2*2));
    phase_shift = [repmat([0,1/2],1,k_left),  zeros(1,N-k_left*2)];    

    rhot = @(t) rho_b./(1 + eps_modulation.*cos(Omega*t + pi*phase_shift)); 

    sqrtkappat = @(t) sqrt(kappa_b)*ones(N,1); % ./[1 + eps*sin(t + pi*(1:N))];
    dkappa = @(t) 0*R;
    d2kappa = @(t) 0*R;
    w3 = @(t) 0*(1:N);

    T = 2*pi/Omega;
    steps = 1000;

    % High contrast parameters \delta
    delta=rho_b/rho0;

    
    % Plot the geometry
    cy = zeros(1,N); cz = zeros(1,N);
    figure, hold on
    t = linspace(0,2*pi);
    for n = 1:N
        plot(cx(n)+R(n)*cos(t), cy(n)+R(n)*sin(t),'k')
    end

    daspect([1 1 1])
    hold off
    close
 
    %% Compute the resonances using the capacitance matrix, which is approximated using the multipole expansion method
    C = capacitance(cx,R,rho0,rho_b,kappa0,kappa_b,delta);

    D = 4/3*pi*R.^3';

    CoeffMat = @(t) makeM(t,delta,kappa0,rho0,D,C,rhot, sqrtkappat, w3);

    %% Solve for Psi

    [Tspan, X_fundamental] = HillSolver(CoeffMat,T,steps); % X_fundamental = [Psi; dPsi/dt]

    %% Solve Floquet exp

    X_T = squeeze(X_fundamental(end,:,:));
    [V,d] = eig(X_T,'vector');
    [f_exp,ind] = sort(real(log(d)/1i/T));
    V = V(:,ind);
%     f_exp_im = imag(log(d(ind))/1i/T);
end