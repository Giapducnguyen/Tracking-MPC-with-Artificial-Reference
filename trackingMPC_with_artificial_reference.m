% % Simulation file for the paper:
% % " MPC for tracking piecewise constant references for constrained
% %   linear systems" - D. Limon et al.

close all;
clear all; 
clc

%% System matrices from the example in Section 4.
A = [1, 1; 0, 1];
B = [0, 0.5; 1, 0.5];
C = [1, 0];

nx = size(A,2);
nu = size(B,2);
ny = size(C,1);
D = zeros(ny,nu);

% % Calculation of the stabilizing gain K and terminal weight P
Q = eye(nx); R = eye(nu);
[K,P,~] = dlqr(A,B,Q,R);
K = -K;
T = 100.*P;

%% Characterization of the steady states and inputs Section 2.1.
A_s = [A-eye(nx),   B, zeros(nx,ny);
        C,           D,      -eye(ny)];
null_A = null(A_s,"rational");
M_theta = null_A(1:(nx+nu),:);
N_theta = null_A(nx+nu+1:end,:);

%% Calculation of invariant set for tracking - Section 2.2.
L = [-K, eye(nu)]*M_theta;
A_w = [A+B*K,           B*L;
       zeros(nu, nx),   eye(nu)];
lambda = 0.99;

% Hard constraints
x_max = 5; u_max = 0.3;
X_set = Polyhedron('lb', -[x_max; x_max], 'ub', [x_max; x_max]);
U_set = Polyhedron('lb', -[u_max; u_max], 'ub', [u_max; u_max]);

z_max = [x_max; x_max; u_max; u_max];
Z_set = Polyhedron('lb', -z_max,'ub',z_max);

%%%%-----------------------%%%%------------------------------%%%%

% % Select 1 of the following computation method

% % Method 1 - Direct computation with tightened constraints:
% O_{\infty,lambda}^w = O(A_w, Wlambda_set)

%{
% % Set constraint on the extended state w = (x, theta)
A_Wlambda = [Z_set.A*[eye(nx), zeros(nx,nu); K, L];
             [zeros(size(Z_set.A,1),nx), Z_set.A*M_theta]];
b_Wlambda = [Z_set.b; lambda.*Z_set.b];
Wlambda_set = Polyhedron('A', A_Wlambda, 'b', b_Wlambda);
Wlambda_set.minHRep();

% % Computation of the maximal admissible invariant set O_{\infty}^w
iter_max = 20; % maximum allowable iterations
found = false; % flag
k = 0; O_setk = Wlambda_set; % initialization

figure(1); hold on; grid on;
title("Computation Evolution of $X_{f} = \mathrm{Proj}_x(\mathcal{O}_{\infty,\lambda}^{w})$","Interpreter","latex");
plot(O_setk.projection(1:nx), 'color', 'black', 'alpha', 0.01, 'edgecolor', 'black');

while (~found && k < iter_max)
    X = sprintf('Iteration: %d',k+1);
    disp(X);

    O_setTemp = Polyhedron('A', Wlambda_set.A * A_w^(k+1), 'b', Wlambda_set.b);
    O_setTemp.minHRep();
    O_setkp1 = intersect(O_setk, O_setTemp);
    O_setkp1.minHRep();
    
    pause(0.5);

    if O_setkp1 == O_setk
        found = true;
        pltf = plot(O_setkp1.projection(1:nx), 'color', 'green', 'alpha', 0.01, 'edgecolor', 'green');
    else
        k = k + 1;
        O_setk = O_setkp1;
        plot(O_setkp1.projection(1:nx), 'color', 'black', 'alpha', 0.01, 'edgecolor', 'black');
    end
end
legend(pltf,"$X_{f} = \mathrm{Proj}_x(\mathcal{O}_{\infty,\lambda}^{w})$","Interpreter","latex");
xlabel("$x_1$","Interpreter","latex");
ylabel("$x_2$","Interpreter","latex");
hold off;
O_setInfty_w = O_setkp1;
%}

% % Method 2 - Indirect computation with transformation and additional constraints

% % following the standard formulation of the paper:
% % "Linear systems with state and control constraints: the theory 
% %  and application of maximal output admissible sets." - E. G. Gilbert

%
% Step 1: System matrix transformation - Equations (2.10 & 4.1) - E. G. Gilbert
U = [eye(nx), (A+B*K - eye(nx))^-1*(-B*L); zeros(nu), eye(nu)];
A_wU = U^-1*A_w*U;
A_wS = A_wU(1:nx, 1:nx);
C_wU = eye(nx+nu)*U;
C_wS = C_wU(:, 1:nx);
C_wL = C_wU(:, nx+1: end);

% Step 2: Set constraint on the extended state (x, theta) - Equation (5.4) - E. G. Gilbert
A_W1 = [Z_set.A*[eye(nx), zeros(nx,nu); K, L];
        Z_set.A*[zeros(size(M_theta,1),nx), M_theta]];
b_W1 = [Z_set.b; Z_set.b];
W1_set = Polyhedron('A', A_W1, 'b', b_W1);
W1_set.minHRep();

A_Wlambda = [Z_set.A*[eye(nx), zeros(nx,nu); K, L];
             Z_set.A*[zeros(size(M_theta,1),nx), M_theta]];
b_Wlambda = [Z_set.b; lambda.*Z_set.b];
Wlambda_set = Polyhedron('A', A_Wlambda, 'b', b_Wlambda);
Wlambda_set.minHRep();

% Step 3: Compute the first set for initialization - Equation (5.4) - E. G. Gilbert
O_set0_A = [W1_set.A*[C_wS, C_wL]; 
            Wlambda_set.A*[zeros(size(C_wL,1),nx), C_wL]]; 
O_set0_b = [W1_set.b; Wlambda_set.b];
O_set0 = Polyhedron('A', O_set0_A, 'b', O_set0_b);
O_set0.minHRep();

figure(1); hold on; grid on;
title("Computation Evolution of $X_{f} = \mathrm{Proj}_x(\mathcal{O}_{\infty,\lambda}^{w})$","Interpreter","latex");
plot(O_set0.projection(1:nx), 'color', 'black', 'alpha', 0.01, 'edgecolor', 'black');

% Step 4: Computation of the maximal admissible invariant set 
% - Equation (5.3) - E. G. Gilbert

% O_{\infty}^w = O_{\infty}( A_wU, \hat{C}_wU, W1_set x Wlambda_set)

% where \hat{C}_wU = [C_wS, C_wL; 0, C_wL]

iter_max = 30;
k = 0; O_setk = O_set0;
found = false;

while(~found && k < iter_max)
    X = sprintf('Iteration: %d',k+1);
    disp(X);

    O_setTemp_A = [W1_set.A*[C_wS*A_wS^(k+1), C_wL]; 
                   Wlambda_set.A*[zeros(size(C_wL,1),nx), C_wL]]; 
    O_setTemp_b = [W1_set.b; Wlambda_set.b];
    O_setTemp = Polyhedron('A', O_setTemp_A, 'b', O_setTemp_b);
    O_setTemp.minHRep();
    O_setkp1 = O_setk.intersect(O_setTemp);
    O_setkp1.minHRep();
    pause(0.5);

    if O_setkp1 == O_setk
        found = true;
    else
        O_setk = O_setkp1;
        k = k + 1;
        plot(O_setkp1.projection(1:nx), 'color', 'black', 'alpha', 0.01, 'edgecolor', 'black');
    end
end

O_setInfty_w = U*O_setkp1; % Equation (2.10) - E. G. Gilbert
O_setInfty_w.minHRep();
pltf = plot(O_setInfty_w.projection(1:nx), 'color', 'green', 'alpha', 0.01, 'edgecolor', 'green');

legend(pltf,"$X_{f} = \mathrm{Proj}_x(\mathcal{O}_{\infty,\lambda}^{w})$","Interpreter","latex");
xlabel("$x_1$","Interpreter","latex");
ylabel("$x_2$","Interpreter","latex");
hold off;
%}

%%%%-----------------------%%%%------------------------------%%%%

%% Calculation of the domain of attraction with prediction horizon N = 3
% % Figure 3 - D. Limon
N = 3;
X_f = O_setInfty_w.projection(1:nx); % Terminal set for state x 
X_setk = X_f; % initialization -> calculation backward

figure(2); hold on; grid on;
title("Computation Evolution of $X_3$","Interpreter","latex");
pltb = plot(X_f, 'color', 'black', 'alpha', 0.01, 'edgecolor', 'black');

for k = 1 : N
    X = sprintf('Iteration: %d',k);
    disp(X);

    set_temp = X_setk + (-B*U_set);
    X_setk = intersect(X_set, Polyhedron('A', set_temp.A*A, 'b', set_temp.b));
    X_setk.minHRep();
    pause(0.5);
    plot(X_setk, 'color', [0.8, 0.8, 0.8], 'alpha', 0.01, 'edgecolor', [0.8, 0.8, 0.8],'Linestyle','--');
end
X_setN = X_setk;
pltf = plot(X_setN, 'color', [0.8, 0.8, 0.8], 'alpha', 0.01, 'edgecolor', 'black','Linestyle','--');
legend([pltb, pltf],["$X_{f}$","$X_3$"],"Interpreter","latex");
xlabel("$x_1$","Interpreter","latex");
ylabel("$x_2$","Interpreter","latex");
hold off;

figure(3);
grid on; hold on;
title("Domain of Attractions and Terminal Sets","Interpreter","latex");
plot(X_setN, 'color', [0.8, 0.8, 0.8], 'alpha', 0.01, 'edgecolor', 'black','Linestyle','--');
plot(X_f, 'color', 'black', 'alpha', 0.01, 'edgecolor', 'black','Linestyle','-');

% % Calculation of O_{\infty}(\hat{x}_s = 0) 
% % and its domain of attraction X_3( O_{\infty}(\hat{x}_s = 0) )
% % Figure 3 - D. Limon

O_set0 = O_setInfty_w.slice([3,4],[0,0]);
X_setk = O_set0;

for k = 1 : N
    X = sprintf('Iteration: %d',k);
    disp(X);

    set_temp = X_setk + (-B*U_set);
    X_setk = intersect(X_set, Polyhedron('A', set_temp.A*A, 'b', set_temp.b));
    X_setk.minHRep();
end
X_setN0 = X_setk;

plot(X_setN0, 'color', 'blue', 'alpha', 0.01, 'edgecolor', 'blue','LineStyle','--');
plot(O_set0, 'color', 'blue', 'alpha', 0.01, 'edgecolor', 'blue','LineStyle','-');
legend({"$X_{3}$","$X_{f}$","$X_3(\mathcal{O}_{\infty}(0))$","$\mathcal{O}_{\infty}(0)$"},"Interpreter","latex");
xlabel("$x_1$","Interpreter","latex");
ylabel("$x_2$","Interpreter","latex");
hold off;

%%%%-----------------------%%%%------------------------------%%%%

%% Closed-loop simulation - Figure 1 & 2 - D. Limon

% % Computed Terminal constraint on (x, theta)
A_Xfw = O_setInfty_w.A;
b_Xfw = O_setInfty_w.b;

% Simulation parameters
Tsim = 90; % [s]
Ts = 1;
Nsim = 90/Ts; % Maximum simulation iterations

% Initialization of decision variable = [state; control; theta]
x0 = [0; -2]; % initial state
x_seq = zeros(nx, N+1); % state sequence
x_seq(:,1) = x0;
u_seq = zeros(nu, N); % control sequence
theta = zeros(nu, 1); % artificial parameter

% Initialization of simulation state, control, target reference, and artificial reference
xSim_seq = nan(nx, Nsim); xSim_seq(:,1) = x0;
uSim_seq = nan(nu, Nsim);
yref_seq = nan(ny, Nsim);
yaRef_seq = nan(ny, Nsim);

% % fmincon options
options = optimoptions("fmincon","Algorithm","interior-point","MaxIterations",100,"Display","iter-detailed");

% Main simulation loop
for k = 1 : Nsim
    % Update current states:
    x0 = xSim_seq(:,k);

    % Update reference
    if k*Ts <= 30
        yt = 4.95;
    elseif k*Ts <= 60
        yt = -5.5;
    else
        yt = 2;
    end
    yref_seq(:,k) = yt;

    % Initialize decision variable for warm-start
    x_seq0 = [reshape(x_seq, nx*(N+1), 1); 
              reshape(u_seq, nu*N, 1);
              theta];

    % Solve for optimal control
    [xi_seq_opt, fval, exitflag, output] = fmincon(...
        @(xi_seq) getCost(xi_seq, yt, N, Q, R, P, T, A, B, M_theta),... % fun
        x_seq0,...                                                      % x0
        [], [], [], [], [], [], ...                                     % A, b, Aeq, beq, lb, ub
        @(xi_seq) getConstraints(xi_seq, x0, N, A, B, x_max, u_max, A_Xfw, b_Xfw),... % nonlcon
        options);                                                       % options

    % Regenerate state sequence from the optimal decision variable
    x_seqtemp = reshape(xi_seq_opt(1:nx*(N+1)), nx, N+1);
    
    % Regenerate control sequence from the optimal decision variable
    u_seq = reshape(xi_seq_opt(nx*(N+1)+1:nx*(N+1)+nu*N), nu, N);
    
    % Regenerate artificial reference from the optimal decision variable
    theta = xi_seq_opt(end-nu+1:end);
    z_s = M_theta*theta;
    x_s = z_s(1:nx);
    u_s = z_s(nx+1:end);
    y_s = N_theta*theta;

    % Store optimal control
    uSim_seq(:,k) = u_seq(:,1);
    % Store artificial reference
    yaRef_seq(:,k) = y_s;

    % Shift solution 1 step for next iteration (warm-start)
    x_seq = [x_seqtemp(:,2:N+1), x_seqtemp(:,N+1)];
    u_seq = [u_seq(:,2:N), K*(x_seqtemp(:,N+1)-x_s) + u_s];

    % Propagate system
    xSim_seq(:,k+1) = A*xSim_seq(:,k) + B*uSim_seq(:,k);

end

%% Plots
t_seq = 0:Ts:Tsim;
figure;
grid on; hold on;
plot(t_seq, xSim_seq(1,:),'-b','LineWidth',1);
plot(t_seq, [nan,yref_seq],'--k','LineWidth',1);
plot(t_seq, [nan,yaRef_seq],'--r','LineWidth',1);
legend({"Actual","Target","Artificial"},"Interpreter","latex");
yticks([-5, 0, 2, 5]);
xlabel("Time","Interpreter","latex");
ylabel("$y$","Interpreter","latex");
title("State Tracking","Interpreter","latex");
hold off;

figure;
grid on; hold on;
pltb = plot(X_f, 'color', [0.8, 0.8, 0.8], 'alpha', 0.2, 'edgecolor', 'black','Linestyle','-.','LineWidth',1);
pltf = plot(X_setN, 'color', [0.8, 0.8, 0.8], 'alpha', 0.01, 'edgecolor', 'black','Linestyle','--','LineWidth',1);
pltx = plot(xSim_seq(1,:), xSim_seq(2,:),'-b','LineWidth',1);
legend([pltb, pltf, pltx],["$X_{f}$","$X_3$", "Actual"],"Interpreter","latex");
xlabel("$x_1$","Interpreter","latex");
ylabel("$x_2$","Interpreter","latex");
title("State Portrait","Interpreter","latex");
hold off;

%%%%-----------------------%%%%------------------------------%%%%

%% Cost function
function V_N = getCost(xi_seq, yt, N, Q, R, P, T, A, B, M_theta)

% Dimensions
nx = size(A,2);
nu = size(B,2);

% Regenerate state sequence from the decision variable
x_seq = reshape(xi_seq(1:nx*(N+1)), nx, N+1);

% Regenerate control sequence from the decision variable
u_seq = reshape(xi_seq(nx*(N+1)+1:nx*(N+1)+nu*N), nu, N);

% Regenerate artificial reference from the decision variable
theta = xi_seq(end-nu+1:end);
z_s = M_theta*theta;
x_s = z_s(1:nx);
u_s = z_s(nx+1:end);
% y_s = N_theta*theta;

% Initialize cost value
V_N = 0;
% Compute cost over prediction horizon
for k = 1 : N
    % Term 1: state tracking
    term1 = (x_seq(:,k)-x_s)'*Q*(x_seq(:,k)-x_s);
    % Term 2: control tracking
    term2 = (u_seq(:,k)-u_s)'*R*(u_seq(:,k)-u_s);
    V_N = V_N + term1 + term2;
end
% Term 3: terminal state tracking
term3 = (x_seq(:,N+1)-x_s)'*P*(x_seq(:,N+1)-x_s);
% Term 4: artificial steady state deviation from target reference
term4 = (x_s-[yt;0])'*T*(x_s-[yt;0]);
% Total cost
V_N = V_N + term3 + term4;
end

%%%%-----------------------%%%%------------------------------%%%%

%% Constraint function
function [c_ineq, c_eq] = getConstraints(xi_seq, x0, N, A, B, x_max, u_max, A_Xfw, b_Xfw)
% Dimensions
nx = size(A,2);
nu = size(B,2);
% Regenerate state sequence from the decision variable
x_seq = reshape(xi_seq(1:nx*(N+1)), nx, N+1);
% Regenerate control sequence from the decision variable
u_seq = reshape(xi_seq(nx*(N+1)+1:nx*(N+1)+nu*N), nu, N);
% Regenerate artificial reference from the decision variable
theta = xi_seq(end-nu+1:end);

% State continuity constraints
c_eq_sc = nan(nx*(N+1),1);
c_eq_sc(1:nx,1) = x_seq(:,1)-x0; 
for k = 1 : N
    xkp1 = A*x_seq(:,k) + B*u_seq(:,k);
    c_eq_sc(nx*k+1:nx*k+nx,1) = x_seq(:,k+1)-xkp1;
end

% State hard constraints from 0 to N-1
lb_x = -[x_max; x_max];
ub_x = [x_max; x_max];
c_ineq_x = nan(2*nx*N,1);
for k = 1 : N
    c_ineq_x(2*nx*(k-1)+1 : 2*nx*(k-1)+2*nx) = [lb_x-x_seq(:,k);
                                                x_seq(:,k)-ub_x];
end

% Control hard constraints from 0 to N-1
lb_u = -[u_max; u_max];
ub_u = [u_max; u_max];
c_ineq_u = nan(2*nu*N, 1);
for k = 1 : N
    c_ineq_u(2*nu*(k-1)+1 : 2*nu*(k-1)+2*nu) = [lb_u-u_seq(:,k);
                                                u_seq(:,k)-ub_u];
end

% Terminal constraint at N: 
% (x(N), theta) \in X_f^w
% X_f^w = {w=(x,theta) | A_Xfw * w <= b_Xfw}
c_ineq_w = A_Xfw * [x_seq(:,N+1);theta] - b_Xfw;

% % Concatenate constraints
% Inequality constraints
c_ineq = [c_ineq_x; % hard constraint on state(from 0 to N-1)
          c_ineq_u; % hard constraint on control(from 0 to N-1)
          c_ineq_w]; % Terminal constraint on state and artificial reference
% Equality constraint
c_eq = c_eq_sc; % Continuity constraints (Multiple shooting method)

end