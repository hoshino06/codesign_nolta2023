% Integrated Design of Physical Structure and Controller 

% Requisite: Symbolic Math Toolbox
% Dependency: fun_galerkin_approximation.m, fun_galerkin_w_sensitivity.m

clear; clc;

% System equations
x = sym('x',[1 2]).';
f = @(rho)[ ( -x(1)^3 - x(2) ); ( x(1)+x(2) ) ];
g = @(rho)[ 0; rho];

% Initial control law derived by feedback linearization
z = [ x(1) ; ( -x(1)^3 +x(2) )];
v_init = 0.4142*z(1) -1.3522*z(2);
u_init = 3*x(1)^5 + 3*x(1)^2*x(2) -x(2) + v_init;

% Parameters of design problem
alpha = 1; 
beta  = 4; 
Js = @(rho) rho;
w     = 1/4; % weighting function
rho   = 1; % initial parameter
u     = u_init; % initial controller
rate  = 0.05; % update rate of rho 

% Basis functions
Phi3 = [x(1)^2, x(1)*x(2), x(2)^2
        ];
Phi8 = [x(1)^2, x(1)*x(2), x(2)^2, ...
        x(1)^4, x(1)^3*x(2), x(1)^2*x(2)^2, x(1)*x(2)^3, x(2)^4
        ];
Phi15 = [ x(1)^2, x(1)*x(2), x(2)^2, ...
          x(1)^4, x(1)^3*x(2), x(1)^2*x(2)^2, x(1)*x(2)^3, x(2)^4, ...
          x(1)^6, x(1)^5*x(2), x(1)^4*x(2)^2, x(1)^3*x(2)^3, ...
          x(1)^2*x(2)^4, x(1)^1*x(2)^5, x(1)^0*x(2)^6
        ];
Phi = Phi15;
N = numel(Phi);

% calculation of the vector W
W = zeros(N,1);
for n = 1:N
    phi_n = Phi(n);
    fun2 = matlabFunction( w * (phi_n), 'Vars',x );
    W(n,1) = integral2( fun2, -1,1, -1,1 );
end

%% Optimization Iteration
num_iter = 1+30;
CN_list  = zeros(num_iter, N);
rho_list = zeros(num_iter, 1);
J_list   = zeros(num_iter, 3);

for i=1:num_iter
    disp(['iteration: ' num2str(i) '/' num2str(num_iter)])

    % Galerkin approximation
    l = x.' * x + u^2;
    [C_N, A, dA, B] = fun_galerkin_w_sensitivity(x,f,g,rho,u,l,Phi);
    CN_list(i,:) = C_N;

    % Current value of objective function
    J_s = Js(rho);
    J_c = (W.')*C_N;
    J  = alpha*J_s + beta*J_c;
    rho_list(i) = rho;
    J_list(i,:) = [J J_s J_c];
    disp(['rho=' num2str(rho) ', J=' num2str(J)])

    % Update of input u
    V_N   = Phi*C_N;
    u = -(1/2)* g(rho)'*jacobian(V_N, x).';

    % Update of design parameter rho
    r = sym('rho');
    dJs = matlabFunction( jacobian(Js(r),r), 'Vars',r);
    dFdrho = alpha* dJs(rho) + beta*(W.')*(A\dA)*(A\B);
    %dFdrho = 0;
    rho = rho - rate*dFdrho;

end

save('results.mat', ...
    "x","f","g","u_init", ... % system setting
    "alpha","beta","w","Js", ... % problem setting 
    "Phi", "N", "num_iter","rate", ... % optimization setting
    "CN_list", "J_list", "rho_list" ... % results
    )

plot(0:num_iter-1, J_list(:,2), '-x' ) 
xlabel('Iteration')
ylabel('Objective function')
