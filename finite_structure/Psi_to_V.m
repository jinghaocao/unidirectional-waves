function V = Psi_to_V(Tspan, Psi, sqrtkappat, D)
%PSI_TO_V transforms Psi (solution of Psi''+ M(t)Psi = 0 to V, given by eq. after eq. (11.4.7) on page 267 in Erik's thesis)
%   INPUT:
%           Tspan       column vector       [t_1, \ldots, t_n] indicating for which times Psi is given
%           Psi         Array               n x N x 2*N array, Psi(k,l,m)
%                                           first coordinate indicates the time t_k,
%                                           the second coordinate l indicates the coordinate l of vector Psi
%                                           the third coordinate m indicates the initial condition e_m for Psi
%                                           Psi(k,l,m) = lth coordinate of Psi at time t_k, when Psi(0) = e_m
%           sqrtkappat  function handle     COLUMN vector valued function handle, sqrt(kappa(t))_i is given by the square root
%                                           of the i^th entry of kappa, which is the (time-dep.) bulk modulus of the ith resonator
%           D           column vector       N-dim vector, D_i indicates the volume of the ith resonator
%   OUTPUT:
%           V           Array               n x N x 2*N array, V(k,l,m)
%                                           first coordinate indicates the time t_k,
%                                           the second coordinate l indicates the coordinate l of vector V i.e V_l in Erik's thesis
%                                           the third coordinate m indicates the initial condition e_m for Psi, 
%                                               that is m indicates the mth basis vector for the solution space of th V vectors
%                                           V(k,l,m) = lth coordinate of V (corresponding to lth resonator) at time t_k, when V(0) = sqrtkappa(0).*D



sqrtD = @(t) (sqrtkappat(t).*D)'; % row vector with entries [sqrt(kappa_1(t))|D_1|, ..., sqrt(kappa_N(t))|D_N|]

N = size(D,1);
steps = size(Tspan,1);

V = NaN(size(Psi));

for m = (1:2*N)
    for k = (1:steps)
        V(k,:,m) = Psi(k,:,m).*sqrtD(Tspan(k));
    end
end

end

