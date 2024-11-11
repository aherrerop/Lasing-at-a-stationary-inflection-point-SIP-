function [x] = SolveLUSystem(A,b)
    % Solve Ax=b by PA = LU decomposition:
    % P*A*x = P*b
    % L*U*x = P*b -> Ux = z, Lz = Pb
    [L,U,P] = lu(A);
    z2 = L \ (P*b);
    x2 = U \ z2;
    
    % Forward Substitution (FS): For Lower triangular matrices:
    [rows,columns] = size(A);
    pb = P*b;
    z=zeros(1,columns);
    z(1,1)=pb(1)./L(1,1);
    %bacward substitution
    for k=2:columns
        z1=(1/L(k,k)).*(pb(k)-sum(L(k,k-1:-1:1).*z(k-1:-1:1)));
        z(1,k)=z1;
    end
    z=conj(z');


    % Backward Substitution (BS)
    x=zeros(1,rows);
    x(1,columns)=z(end)./U(columns,columns);    
    %bacward substitution
    for k=columns-1:-1:1
        x1=1/U(k,k).*(z(k)-sum(U(k,k+1:end).*x(k+1:end)));
        x(1,k)=x1;
    end
    x=conj(x');
end