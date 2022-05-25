function [TOUT, X_fund] = HillSolver(CoeffMat,T,steps)
    %HillSolver
    %   Output: returns array: fundamental solution of degree one equation associated to (11.4.10) in Erik's thesis for one period i.e. on [0,T]
    %   Input: 
    %       CoeffMat    function handle     coefficient matrix CoeffMat = M of Hill
    %                                       equation Psi'' + M(t)Psi = 0 ((11.4.10) in Erik's thesis), which is T periodic
    %       T           scalar              period of CoeffMat
    %       steps       integer > 0         number of steps where the solution
    %                                       Psi will be evaluated (Psii will be evalutated on linspace(0,T,steps))
    %   Output:
    %       Tspan       Column vector       Column vector indicating where Psi
    %                                       has been evaluate (Tspan = linspace(0,T,steps)')
    %       X_fund      Array               steps x 2*N x 2*N, where CoeffMat(t) is a
    %                                       N x N matrix. 
    %                                       X_fund(m,n,k) = (Psi(Tspan(m)))_n, where
    %                                       Psi is solution to the above Hill equation with initial condition e_k
    %                                       X_fund(m,:,:) is the 1 x 2*N x 2*N solution
    %                                       of the (first order) Hill equation at time t_m with initial cond. Id_(2*N)  
    
 
    N = size(CoeffMat(0),2);   
    Id = eye(N);
    Z = zeros(N);

    odefun = @(t,X) [Z, Id; -CoeffMat(t), Z]*X;

    Tspan = linspace(0,T,steps)';

    X_fund = NaN(steps,2*N,2*N); % X_fund(m,n,l) = (X_l)_n(t(m)) i.e. nth entry of solution X at time t(m) of ODE with initial condition e_l
    
    Id_2N = eye(2*N);

    for n = (1:2*N)
        [TOUT,X_fund(:,:,n)] = ode45(odefun,Tspan,Id_2N(:,n));
    end
    
end