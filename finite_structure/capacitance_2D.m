function C = capacitance_2D(c,R,rho0,rho_b,kappa0,kappa_b,delta)
cx = c(:,1);
cy = c(:,2);
N = length(cx);
%     N_multi = 4;
N_multi = 0;

z = delta;
A = MakeA(R,z,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy);

B = zeros(N);

for i = 1:N
    B(i,i) = A(2*i-1,2*i-1);
end

for i = 1:N
    for j = 1:N
        if i ~= j
            B(i,j) = -A(2*i-1,2*j);
        end
    end
end

B = -B/4/pi/R(1)^2; % surface area of a sphere

C = B^-1;