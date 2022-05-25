N = 13;
eps = 0.2;

[V_mode,f_exp,ind] = produce_fexp(N,eps);

figure
hold on
plot(1:N*2,f_exp(1:N*2),'.k')


% find out edge modes ind
ind_1 = (N+1)/2;
ind_2 = ind_1+N;
figure
hold on
for i = 1:N*2
    subplot(5,ceil(N*2/5),i)
    plot(1:N,real(V_mode(1:N,i)))
    xlim([1,N])
end
figure
plot(1:N,real(V_mode(1:N,ind(ind_1))))
xlim([1,N])
figure
plot(1:N,real(V_mode(1:N,ind(ind_2))))
xlim([1,N])

function [V_mode,f_exp,ind] = produce_fexp(N,eps)
    % Number of Resonators
    % radii
    R = 0.1*ones(1,N);
    Omega = 0.3;
    %%% Material parameters
    rho0 = 1e3;             % background material
    kappa0 = 1e3;           % background material
    v = sqrt(kappa0/rho0);  % speed of sound in water

    rho_b = 1;            % density of resonators  
    kappa_b = 1;          % bulk modulus of resonators
    v_b = sqrt(kappa_b/rho_b);  % speed of sound resonators
    
    d_1 = 1.5;
    u_1 = [0 0.5];
    cx = u_1;
    for j = 1:floor(N/2)/2-1
        cx = [cx u_1+d_1*j];
    end
    cx = [cx cx(end)+1];
    middle = cx(end);
    for j = 1:floor(N/2)/2
        cx = [cx u_1+middle+j*d_1-0.5];
    end
    
%     cx = 0.5*[0:1:N-1];
%     cx = cx + [zeros(1,ceil(N/2)) 0.25*ones(1,N-ceil(N/2))];
    % Plot the geometry
    cy = [zeros(1,N)]; 
    cz = zeros(1,N);
    figure, hold on
    t = linspace(0,2*pi);

    for n = 1:N
        plot(cx(n)+R(n)*cos(t), cy(n)+R(n)*sin(t),'k')
    end
%     cx = flip(cx);
    daspect([1 1 1])
    hold off
%     close
    c = [cx;cy]';
    eps_modulation =[];
    phase_shift = [];
    N_modulated = floor(N/2);
    eps_modulation = [ones(1,N_modulated)*eps zeros(1,N-N_modulated)];
    phase_shift = [repmat([1/3,1/2],1,N_modulated/2) zeros(1,N-N_modulated)];
    rhot = @(t) rho_b./(1 + eps_modulation.*cos(Omega*t + pi*phase_shift)); 

    sqrtkappat = @(t) sqrt(kappa_b)*ones(N,1); % ./[1 + eps*sin(t + pi*(1:N))];
    dkappa = @(t) 0*R;
    d2kappa = @(t) 0*R;
    w3 = @(t) 0*(1:N);

    T = 2*pi/Omega;
    steps = 100;

    % High contrast parameters \delta
    delta=rho_b/rho0;

    

    %% Compute the resonances using the capacitance matrix, which is approximated using the multipole expansion method
    C = capacitance(cx,R,rho0,rho_b,kappa0,kappa_b,delta);

    D = 4/3*pi*R';

    CoeffMat = @(t) makeM(t,delta,kappa0,rho0,D,C,rhot, sqrtkappat, w3);

    %% Solve for Psi

    [Tspan, X_fundamental] = HillSolver(CoeffMat,T,steps); % X_fundamental = [Psi; dPsi/dt]
    [TOUT, X_bloch,V_mode,V,D] = Hill2BlochModes(CoeffMat,T,steps);
    
    %% Solve Floquet exp
    
    X_T = squeeze(X_fundamental(end,:,:));
    [V,d] = eig(X_T,'vector');
    [f_exp,ind] = sort(real(log(d)/1i/T));
    f_exp_im = imag(log(d(ind))/1i/T);

end