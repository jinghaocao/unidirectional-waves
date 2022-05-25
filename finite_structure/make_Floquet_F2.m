% The folding numbers are zero
function F2 = make_Floquet_F2(A0,A1_f,A2_f,F1,T,fn,p_index)
   n = length(F1);
   F2 = zeros(n);
   for k = p_index:p_index+1
       for l = p_index:p_index+1
            for j = 1:n
                for m = [-1 0 1]
                    if m == fn(j)-fn(l)
                        0;
                    else
                        if abs(fn(k)-fn(l)-m)>1
                            a_1 = 0;
                        else
                            a_1 = A1_f(fn(k)-fn(l)-m+2,k,j);
                        end
                        if abs(fn(k)-fn(j))>1
                            a_2 = 0;
                        else
                            a_2 = A1_f(fn(k)-fn(j)+2,k,j);
                        end
                    F2(k,l) = F2(k,l) + (a_1-a_2)*A1_f(m+2,j,l)/(2*pi*1i*m/T +A0(l,l)-A0(j,j));
                    k
                    l
                    j
                    end
               end
            end
            for j = 1:n
                for m = [-1 0 1]
                    if m == fn(k)-fn(j)
                        0;
                    else
                        F2(k,l) = F2(k,l) + A1_f(m+2,k,j)*F1(j,l)/((2*pi*1i*m/T +A0(j,j)-A0(k,k)));
                    end
                end
            end
            if abs(fn(k)-fn(l))<4
                F2(k,l)= F2(k,l)+A2_f(k,l,fn(k)-fn(l)+3);
            end
%             if F2(k,l)<0.0000001
%                 F2(k,l) =0;
%             end
       end
   end
end