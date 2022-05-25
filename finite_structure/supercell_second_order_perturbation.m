function [F2,v_1_1,v_1_2,beta,rate]=supercell_second_order_perturbation(M_0,edge_ind,Omega,N,phase_shift)    
% Solve perturbation analytically
%     M_0 = CoeffMat(0);
%     M_0 = delta*kappa0/rho0*inv(diag(D))*C;
    A0 = [zeros(N) eye(N); -M_0 zeros(N)];
    T = 2*pi/Omega;    
% diagonalize A0 and keep the basis
    [Q,d] = eig(A0,'vector');
    % compute the folding number fn and set up F_0
    fn = round(imag(d)/Omega);
    F0 = diag(real(d)+1i*(imag(d)-(round(imag(d)/Omega))*Omega));
    f_0 = diag(real(F0/(1i)));
    [f_0,index]=sort(f_0,'descend');
    f_0(edge_ind(2))=f_0(edge_ind(1));
    F0 = 1i*diag(f_0);
    fn = fn(index);
    % sort basis transformation matrix
    Q = Q(:,index);
    A0 = Q\A0*Q;
    A0(edge_ind(2),edge_ind(2))=A0(edge_ind(1),edge_ind(1));
    % constant fourier matrix of A1 is zero
    A1_f = zeros(3,2*N,2*N);
    % 1st and -1st fourier matrix of A1
    M = M_0;
    M1 = zeros(N);
    M2 = zeros(N);
    for l = 1 : N
        for j = 1: N
            M1(l,j)=0.5*M(l,j)*(exp(1i*phase_shift(j)-exp(1i*phase_shift(l))));
            M2(l,j)=0.5*M(l,j)*(exp(-1i*phase_shift(j))-exp(-1i*phase_shift(l)));
        end
    end
    A1_f(1,:,:) = Q\[zeros(N) zeros(N); -M2 zeros(N)]*Q;
    A1_f(3,:,:) = Q\[zeros(N) zeros(N); -M1 zeros(N)]*Q;
    F1 = make_Floquet_F1(F0,A0,A1_f,T,fn);
%   project index
%   p_index = 2*N-ind_deg_1;
    p_index = edge_ind(1);
    % Projection matrix P
    P = zeros(2*N,2*N);
    P(p_index,p_index)=1;
    P(p_index+1,p_index+1)=1;
    % Make A2_0
    M2_0 = zeros(N);
    M2_2 = zeros(N);
    M2_m2 = zeros(N);
    for l = 1 : N
        for j = 1: N
            M2_0(l,j)=M(l,j)*(0.5-0.25*(exp(1i*(phase_shift(j)-phase_shift(l)))+exp(1i*(phase_shift(l)-phase_shift(j)))));
            M2_2(l,j)=M(l,j)*(-0.25*exp(1i*(phase_shift(l)+phase_shift(j)))+0.25*exp(1i*2*phase_shift(l)));
            M2_m2(l,j)=M(l,j)*(-0.25*exp(-1i*(phase_shift(l)+phase_shift(j)))+0.25*exp(-1i*2*phase_shift(l)));
        end    
    end
    A2_f = zeros(2*N,2*N,5);
    A2_f(:,:,3) = Q\[zeros(N) zeros(N); -M2_0 zeros(N)]*Q;
    A2_f(:,:,1) = Q\[zeros(N) zeros(N); -M2_m2 zeros(N)]*Q;
    A2_f(:,:,5) = Q\[zeros(N) zeros(N); -M2_2 zeros(N)]*Q;

    F_2 = make_Floquet_F2(A0,A1_f,A2_f,F1,T,fn,p_index);
    F2 =  F_2(p_index:p_index+1,p_index:p_index+1)
    % effective hamiltonian
    G = zeros(2*N);
    for k = [1:p_index-1,p_index+2:2*N]
        G(k,k)=(f_0(p_index)-f_0(k))^(-1);
    end
    H = P*F1*G*F1*P;
    H = H(p_index:p_index+1,p_index:p_index+1)+F2;
    [beta,rate] = eig(H);
    v_0_1 = zeros(2*N,1);
    v_0_1(p_index:p_index+1)=beta(:,1);
    v_0_2 = zeros(2*N,1);
    v_0_2(p_index:p_index+1)=beta(:,2);
    
    v_1_1 = Q*G*F1*v_0_1;
    v_1_2 = Q*G*F1*v_0_2;


% Thea's implementation gives the same results    
%     % constant fourier matrix of A1 is zero
%     A1_thea = zeros(2*N,2*N,3);
%     % 1st and -1st fourier matrix of A1
%     M = M_0;
%     M1 = zeros(N);
%     M2 = zeros(N);
%     for l = 1 : N
%         for j = 1: N
%             M1(l,j)=0.5*M(l,j)*(exp(1i*phase_shift(j)-exp(1i*phase_shift(l))));
%             M2(l,j)=0.5*M(l,j)*(exp(-1i*phase_shift(j))-exp(-1i*phase_shift(l)));
%         end
%     end
%     A1_thea(:,:,1) = Q\[zeros(N) zeros(N); -M2 zeros(N)]*Q;
%     A1_thea(:,:,3) = Q\[zeros(N) zeros(N); -M1 zeros(N)]*Q;
%     F2_matrix = thea_make_F2(A0,A1_thea,A2_f,Omega);
%     F2_matrix_thea=F2_matrix(p_index:p_index+1,p_index:p_index+1)
%     rate_thea = eig(F2_matrix_thea);
end    
