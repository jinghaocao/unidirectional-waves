function [TOUT, X_bloch,V_mode,V,D] = Hill2BlochModes(CoeffMat,T,steps)
%HILL2BLOCHMODES
    %   Output: returns array: Bloch modes of degree one equation associated to (11.4.10) in Erik's thesis for one period i.e. on [0,T]
    %   Input: 
    %       CoeffMat    function handle     coefficient matrix CoeffMat = M of Hill
    %                                       equation Psi'' + M(t)Psi = 0 ((11.4.10) in Erik's thesis), which is T periodic
    %       T           scalar              period of CoeffMat
    %       steps       integer > 0         number of steps where the solution
    %                                       Psi will be evaluated (Psii will be evalutated on linspace(0,T,steps))
    %   Output:
    %       Tspan       Column vector       Column vector indicating where Psi
    %                                       has been evaluate (Tspan = linspace(0,T,steps)')
    %       X_bloch      Array               steps x 2*N x 2*N, where CoeffMat(t) is a
    %                                       N x N matrix. 
    %                                       X_bloch(m,n,k) = (Psi(Tspan(m)))_n, where
    %                                       Psi is solution to the above Hill equation with initial condition e_k
    %                                       X_bloch(m,:,:) is the 1 x 2*N x 2*N solution
    %                                       of the (first order) Hill equation at time t_m with initial cond. Id_(2*N)  
    %
    %
    %   Procedure:
    %       1)  Compute fundamental solution X_fund with HillSolver
    %       2)  Diagonalize fundamental solution X_fund -> V = [v_1, \ldots, v_2N] eigenvectors
    %       3)  Compute Bloch modes by solving the initial value problem
    %               Psi'' = [0, Id; -M(t), 0]Psi
    %               Psi(0) = V
    %           with ode45.

    N = size(CoeffMat(0),2);   

    % 1)  Compute fundamental solution X_fund with HillSolver
    [TOUT, X_fund] = HillSolver(CoeffMat,T,steps);
    

    %%%%   CHECK IF THIS IS REALLY CORRECT!!!!!
    % 2)  Diagonalize fundamental solution X_fund -> V = [v_1, \ldots, v_2N] eigenvectors
    [V,D] = eig(reshape(X_fund(end,:,:),[2*N,2*N])); %%%%   CHECK IF THIS IS REALLY CORRECT!!!!!
    %%%%   CHECK IF THIS IS REALLY CORRECT!!!!!


    % 3)  Compute Bloch modes by solving the initial value problem

    Id = eye(N);
    Z = zeros(N);

    odefun = @(t,X) [Z, Id; -CoeffMat(t), Z]*X;

    Tspan = linspace(0,T,steps)';
%     Tspan = linspace(0,2*T,steps)';

    X_bloch = NaN(steps,2*N,2*N); % X_bloch(m,n,l) = (X_l)_n(t(m)) i.e. nth entry of solution X at time t(m) of ODE with initial condition v_l


    for n = (1:2*N)
        [TOUT,X_bloch(:,:,n)] = ode45(odefun,Tspan,V(:,n));
    end

    %%% TEST : whether we really get Bloch modes %%%%%

    V_mode = reshape(X_bloch(end,:,:),[2*N,2*N]);

    if norm(V_mode - V*D) < 0.0001
        if norm(V_mode - reshape(X_bloch(1,:,:),[2*N,2*N])*D) < 0.0001
        disp('all good')
        else
            disp('X_bloch(T) - X_bloch(0)*D too big, approximation is not close enough to actual Bloch modes')
            error_1 = norm(V_mode - reshape(X_bloch(1,:,:),[2*N,2*N])*D);
        end
    else
        disp('X_bloch(T) - eig_vectors(X_fund(T))*D too big, approximation is not close enough to actual Bloch modes')
        error_2 = norm(V_mode - V*D)
    end

end

