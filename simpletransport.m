function simpletransport

% Simulates the transport through nanowire. Note that everything will be
% starting at l = 1, there is no 0th element in matlab, this does not
% matter as everything will simply be shifted by +1.

t=1 ; % value of "t" does not matter, it will cancel itself out, why ????
size_M=2;
size_L=6; 
lambda_f=10;
B=0;
Vl(1:size_L)=zeros(1,6);

%initialize Cl_1_2 vector of matrices
Cl_1 = eye(size_M);
Cl_2 = zeros(size_M);
%initialize channelcoefficients R(+-) and L(+-)
L_plus  = eye(size_M);
L_minus = zeros(size_M);
l=1;
Ef_scalar = calculate_Ef(lambda_f,t);
[Hl_lmin1,Hl_lplus1]=build_Hl_lmin1_lplus1(B,t,size_M);
Tell=Tl(l,Vl,t,Ef_scalar,size_M,Hl_lmin1,Hl_lplus1)
%[U_plus, U_minus, lambda_plus, lambda_minus] = get_U_lambda(Tl,size_M)
%T_zero=build_T_zero(U_plus, U_minus, lambda_plus, lambda_minus,size_M)
[T_zero, lambda] = build_T_zero(Tell,size_M)
%Pl = build_Pl(Tl,Cl_1, Cl_2,size_M);

%T = T_zero\Tl*Tl*T_zero*[L_plus; L_minus]
TN = 1;

 for(j= 1:size_L)
     TN=Tl(l,Vl,t,Ef_scalar,size_M,Hl_lmin1,Hl_lplus1)*TN
     l= l+1;
     
 end
 
 T = T_zero/TN*T_zero
end

function Hl=build_Hl(l,Vl,t,size_M)
    Hl=zeros(size_M,size_M);
    for(i=1:size_M)
        for(m=1:size_M)
            if(i==m)
                Hl(i,m)=Vl(l)+4*t;
            elseif(i==m+1 || i==m-1)
                Hl(i,m)=-t;
            end
        end
    end
end


function [Hl_lmin1,Hl_lplus1]=build_Hl_lmin1_lplus1(B,t,size_M)
    Hl_lmin1=zeros(size_M,size_M);
    for(i=1:size_M)
        for(m=1:size_M)
            if(i==m)
                Hl_lmin1(i,m)=-t*exp(2*pi*1i*B*m);
            end
        end
    end
    Hl_lplus1=conj(Hl_lmin1);
    
end

function [Tl]=Tl(l,Vl,t,Ef_scalar,size_M,Hl_lmin1,Hl_lplus1)
    %Set up matrices needed later
    Ef=eye(size_M,size_M)*Ef_scalar;
    Hl=build_Hl(l,Vl,t,size_M);
    %Initialize Tl and fill with physics
    Tl                                      =zeros(2*size_M,2*size_M);
    Tl(1:size_M,size_M+1:2*size_M)          =eye(size_M);
    Tl(size_M+1:2*size_M,1:size_M)          =-Hl_lplus1\Hl_lmin1;
    Tl(size_M+1:2*size_M,size_M+1:2*size_M) =Hl_lplus1\(Ef-Hl);
end

function Ef_scalar=calculate_Ef(lambda_f, t)
    Ef_scalar=2*t*(1-cos(2*pi/lambda_f));
end

function [T_zero,lambda] = build_T_zero(Tl,size_M)
    [T_zero, lambda] = eig(Tl);
    %normalize in respect to the upper half
    for(i= 1:2*size_M)
        z(i) = sqrt(T_zero(1:size_M,i)'*T_zero(1:size_M,i));
    end
    %z = sqrt(sum(T_zero(1:size_M,:).^2,1));
    T_zero = T_zero./z(ones(size(T_zero,1),1),:);
end

function [Pl] = build_Pl(Tl,Cl_1, Cl_2,size_M)
    Pl      = eye(2*size_M,2*size_M);
    Pl_2    = inv(Tl(size_M+1:2*size_M,1:size_M)*Cl_2+Tl(size_M+1:2*size_M,size_M+1:2*size_M));
    Pl_1    = -Pl_2*Tl(size_M+1:2*size_M,1:size_M)*Cl_1;
    Pl(size_M+1:2*size_M,1:size_M)     = Pl_1;
    Pl(size_M+1:2*size_M,size_M+1:2*size_M) = Pl_2;
end

%function [U_plus, U_minus, lambda_plus, lambda_minus] = get_U_lambda(Tl,size_M)
    %solves eigenvalue equation to generate eigenvectors u_m(+-), + is assumed to be the
    %first half of eigenvectors and - the second. The eigenvalues are
    %labeled analogously.
    %[V, D] = eig(Tl);
    %if(size(V) == 2*size_M)
    %    U_plus  = V(1:size_M,1:size_M);
    %    %U_minus = V(size_M+1:2*size_M,size_M+1:2*size_M);
    %    U_minus = V(1:size_M,size_M+1:2*size_M);
    %    lambda_plus = D(1:size_M,1:size_M);
    %    lambda_minus = D(size_M+1:2*size_M,size_M+1:2*size_M);
    %else
    %   disp('Eigenvalue Matrix dimension mismatch, nequal size_M');
    %end
%end

%function [T_zero] = build_T_zero(U_plus, U_minus, lambda_plus, lambda_minus,size_M)
     %builds T_zero from matrices U(+-) and lambda(+-)
     %T_zero(1:size_M,1:size_M) = U_plus;
     %T_zero(1:size_M,size_M+1:2*size_M) = U_minus;
     %T_zero(size_M+1:2*size_M,1:size_M) = U_plus*lambda_plus;
     %T_zero(size_M+1:2*size_M,size_M+1:2*size_M) = U_minus*lambda_minus;
%end