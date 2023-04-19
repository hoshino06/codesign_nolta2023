function [C_N, A, dAdrho,B] = fun_galerkin_w_sensitivity(x, f, g, rho, u, l, Phi)

% default argument for debug
switch nargin
    case 0
        % System equations
        x = sym('x',[1 2]).';
        f = @(rho)[ ( -x(1)^3 - x(2) ); ( x(1)+x(2) ) ];
        g = @(rho)[ 0; rho];
        % Initial paramter
        rho = 1;
        % Initial Control input
        z = [ x(1) ; ( -x(1)^3 +x(2) )];
        v = 0.4142*z(1) -1.3522*z(2);
        u = 3*x(1)^5 + 3*x(1)^2*x(2) -x(2) + v;
        % Cost function
        l = x.' * x + u^2;
        % Basis functions
        Phi = [ x(1)^2, x(1)*x(2), x(2)^2, ...
                %x(1)^4, x(1)^3*x(2), x(1)^2*x(2)^2, x(1)*x(2)^3, x(2)^4 ...
        ]; 
end

N = numel(Phi);
A = zeros(N);
B = zeros(N,1);
dAdrho = zeros(N);

for n = 1:N
    phi_n = Phi(n);
    for j = 1:N
        phi_j = Phi(j);

        % calculation of the matrix A 
        funA = matlabFunction( ( jacobian(phi_j, x) * (f(rho)+g(rho)*u) ) * (phi_n) );
        A(n,j) = integral2( funA, -1,1, -1,1);

        % calculation of dA/drho
        r = sym('rho');
        df = jacobian(f(r),r);
        dg = jacobian(g(r),r);
        funDA = matlabFunction( ( jacobian(phi_j, x) * ( df+dg*u ) ) * (phi_n), 'Vars', x); 
        if abs(funDA(0,0))+abs(funDA(1,0))+abs(funDA(0,1))+abs(funDA(1,1)) > 0 
        % integral2 fails when fun3 is equivalently eqaul to 0
            dAdrho(n,j) = integral2( funDA, -1,1, -1,1);
        end
    end

    % calculation of b 
    funB = matlabFunction( l * (phi_n) );
    B(n,1) = integral2( funB, -1,1, -1,1 );
end

C_N = - A\B;
