% A1: n*N*N tensor where i*N*N entries are the (i-2)th Fourier coefficients
% fn: array of the folding number
% p stores the labeling of eigenvalues, 0 stands for non-degenerate, 1 in
% position i means that the 1st degenerate eigenvalues are in positions i,
% 2 in j etc...
function F1 = make_Floquet_F1(F0,A0,A1_f,T,fn)
    N = length(F0);
    F1 = zeros(N);
    % Formula for diagonals
    for k = 1:N
        F1(k,k)=A1_f(2,k,k);
    end
    % Formula for off-diagonals
    for k = 1:N
        for l = 1:N
            if abs(F0(k,k)-F0(l,l))<0.0001
               if abs(fn(k)-fn(l))>1
                   F1(k,l) =0;
               else
               F1(k,l) = A1_f(fn(k)-fn(l)+2,k,l);
               end
            else
                sum = 0;
                for m = -1:1
                    sum = sum + A1_f(m+2,k,l)/(2*pi*1i*m/T+A0(l,l)-A0(k,k));
                end
                F1(k,l)=(F0(l,l)-F0(k,k))*sum;
            end
%             if F1(k,l)<0.0000001
%                 F1(k,l) =0;
%             end
        end
    end
end