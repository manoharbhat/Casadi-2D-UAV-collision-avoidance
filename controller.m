classdef controller < matlab.System & matlab.system.mixin.Propagates
% untitled Add summary here

% This template includes the minimum set of functions required
% to define a System object with discrete state.

properties
% Public, tunable properties.

end

properties (DiscreteState)
end

properties (Access = private)
% Pre-computed constants.
casadi_solver 
x0%initial point
lbx%lower limit
ubx%upper limit
lbg
ubg
p%parameters- includes initial and final reference
end

methods (Access = protected)
function num = getNumInputsImpl(~)
    num = 2;
end
function num = getNumOutputsImpl(~)
    num = 1;
end
function dt1 = getOutputDataTypeImpl(~)
    dt1 = 'double';
end
function dt1 = getInputDataTypeImpl(~)
    dt1 = 'double';
end
function sz1 = getOutputSizeImpl(~)
    sz1 = [2,1];
end
function sz1 = getInputSizeImpl(~)
    sz1 = [6,1];
end
function cp1 = isInputComplexImpl(~)
    cp1 = false;
end
function cp1 = isOutputComplexImpl(~)
    cp1 = false;
end
function fz1 = isInputFixedSizeImpl(~)
    fz1 = true;
end
function fz1 = isOutputFixedSizeImpl(~)
    fz1 = true;
end
function setupImpl(obj,~,~)
% Implement tasks that need to be performed only once such as pre-computed constants.

addpath('casadi')
import casadi.*

T = 0.1; % sampling time (s)
N = 10; % prediction horizon

v_max = 10; v_min = -v_max;
omega_max = pi; omega_min = -omega_max;

x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
states = [x;y;theta]; n_states = length(states);

v = SX.sym('v'); omega = SX.sym('omega');%control variables
controls = [v;omega]; n_controls = length(controls);
rhs = [v*cos(theta);v*sin(theta);omega]; % system r.h.s or model equations 
U = SX.sym('U',n_controls,N); % Decision variables (controls)
X = SX.sym('X',n_states,(N+1)); % Decision variable (states)
% A vector that represents the states over the optimization problem.
P = SX.sym('P',n_states + n_states);
% parameters (which include the initial state and the reference state)

f = casadi.Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
 
%states and controls are inputs and rhs is output
J = 0; % cost function
g = [];  % constraints vector

Q = zeros(3,3); Q(1,1) = 1000; Q(2,2) = 1000; Q(3,3) = 100; % weighing matrices (states)
R = zeros(2,2); R(1,1) = 10; R(2,2) = 10; % weighing matrices (controls)

st  = X(:,1); % initial state
g = [g;st-P(1:3)]; % initial condition constraints
for k = 1:N%taking from 0 to N-1 or 1:N
    st = X(:,k);%  states
    con = U(:,k);% controls
    xst = X(1,k);
    yst = X(2,k);
    J = J + (st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con; % calculate cost
    st_next = X(:,k+1);
    f_value = f(st,con);
    st_next_euler = st + (T*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
end
% make the decision variable one column  vector
OPT_variables = [reshape(X,3*(N+1),1);reshape(U,2*N,1)];

nlp_prob = struct('f', J, 'x', OPT_variables, 'g', g, 'p', P);%nlp programming

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);
obj.casadi_solver = solver;

%obj = struct;

obj.lbg(1:3*(N+1)) = 0;  % -1e-20  % Equality constraints
obj.ubg(1:3*(N+1)) = 0;  % 1e-20   % Equality constraints

obj.lbx(1:3:3*(N+1),1) = -12; %state x lower bound
obj.ubx(1:3:3*(N+1),1) = 12; %state x upper bound
obj.lbx(2:3:3*(N+1),1) = -12; %state y lower bound
obj.ubx(2:3:3*(N+1),1) = 12; %state y upper bound
obj.lbx(3:3:3*(N+1),1) = -pi; %state theta lower bound
obj.ubx(3:3:3*(N+1),1) = pi; %state theta upper bound

obj.lbx(3*(N+1)+1:2:3*(N+1)+2*N,1) = v_min; %v lower bound
obj.ubx(3*(N+1)+1:2:3*(N+1)+2*N,1) = v_max; %v upper bound
obj.lbx(3*(N+1)+2:2:3*(N+1)+2*N,1) = omega_min; %omega lower bound
obj.ubx(3*(N+1)+2:2:3*(N+1)+2*N,1) = omega_max; %omega upper bound

end

function u = stepImpl(obj,x,t)
disp(t);
N = 10;   
tic;
xm = x(1:3); % initial condition.
xr = x(4:6); % Reference posture.
u0 = zeros(N,2);        % two control inputs for each robot
X0 = repmat(xm,1,N+1)'; % initialization of the states decision variables
% initial value of the optimization variables
obj.x0  = [reshape(X0',3*(N+1),1);reshape(u0',2*N,1)];
obj.p   = [xm;xr]; % set the values of the parameters vector
solver = obj.casadi_solver;
sol = solver('x0', obj.x0, 'lbx', obj.lbx, 'ubx', obj.ubx,...
        'lbg', obj.lbg, 'ubg', obj.ubg,'p', obj.p);
u = reshape(full(sol.x(3*(N+1)+1:end))',2,N)';
u = u(1,:)';
toc
end
function resetImpl(obj)
    % Initialize discrete-state properties.
end
end
end