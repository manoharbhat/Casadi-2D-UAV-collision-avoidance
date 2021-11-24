classdef momul< matlab.System & matlab.system.mixin.Propagates
    % untitled Add summary here
    %
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
        x0
        lbx
        ubx
        lbg
        ubg
        p
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
        	sz1 = [12,1];
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
        
        %%
        function setupImpl(obj,~,~)  % you used J instead of obj.
            % Implement tasks that need to be performed only once, 
            % such as pre-computed constants.
            
           addpath('casadi')
           import casadi.*
           f_g=1.7715;
t_g=1.1699;


l_a=0.6604;
l_w=0.4699;
l_d=0.05;
l_e=0.02;
l_h=0.1778;

m_h=1.15;
m_w=1.87;
m_f=0.5750;
m_b=0.5750;

k_f=0.1188;



g=9.81;


p_a=atan((l_d+l_e)/l_a);
p_w=atan(l_d/l_w);

j_e=m_w*l_w*l_w+(m_f+m_b)*l_a*l_a;
j_p=(m_h+m_b)*l_h*l_h;
j_l=m_w*l_w*l_w+(m_f+m_b)*(l_h*l_h+l_a*l_a);


   
           
            T = 1; % Time horizon
            N = 10; % number of control intervals
            
           
            v_f_max = 24; 
            v_f_min = -24;
            v_b_max = 24; 
            v_b_min = -24;
 
            % Declare model variables
            x1 = SX.sym('x1');
            x2 = SX.sym('x2');
            x3 = SX.sym('x3');
            x4 = SX.sym('x4'); 
            x5 = SX.sym('x5'); 
            x6 = SX.sym('x6');
            states = [x1;x2;x3;x4;x5;x6];
            n_states = length(states);
            
            v_f = SX.sym('v_f');
            v_b = SX.sym('v_b');
            controls = [v_f;v_b];
            n_controls = length(controls);
            
            %Modell
           
            rhs =[x4; x5; x6; -(k_f*l_a*(v_f+v_b)-l_a*f_g)/j_e; (k_f*l_h*(v_f-v_b))/j_p; (-f_g*sin(x2)*l_a)/j_l];        

U = SX.sym('U',n_controls,N); % Decision variables (controls)
X = SX.sym('X',n_states,(N+1)); % Decision variable (states)
% A vector that represents the states over the optimization problem.
P = SX.sym('P',n_states + n_states);
% parameters (which include the initial state and the reference state)
         

f = casadi.Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)

J = 0; % Objective function
g = [];  % constraints vector


Q = zeros(6,6); Q(1,1) = 500;Q(2,2) = 1;Q(3,3) = 50;Q(4,4) = 1;Q(5,5) = 1 ;Q(6,6) = 10; % weighing matrices (states)
R = zeros(2,2); R(1,1) = 5; R(2,2) = 5; % weighing matrices (control

st  = X(:,1); % initial state
g = [g;st-P(1:6)]; % initial condition constraints

           
     for k = 1:N
    st = X(:,k); 
    con = U(:,k);
    J = J+(st-P(7:12))'*Q*(st-P(7:12)) + con'*R*con; % calculate obj
    st_next = X(:,k+1);
    f_value = f(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
    end
          
% QN = 10*eye(6);
%   J = J+(st-P(7:12))'*Q*(st-P(7:12)) + con'*R*con;
% make the decision variable one column  vector

OPT_variables = [reshape(X,6*(N+1),1);reshape(U,2*N,1)];

nlp_prob = struct('f', J, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

obj.casadi_solver = solver;
       

obj.lbg(1:6*(N+1)) = 0;  % -1e-20  % Equality constraints
obj.ubg(1:6*(N+1)) = 0;  % 1e-20   % Equality constraints

obj.lbx(1:6:6*(N+1),1) = -inf; %state x lower bound
obj.ubx(1:6:6*(N+1),1) = inf; %state x upper bound
obj.lbx(2:6:6*(N+1),1) = -inf; %state y lower bound
obj.ubx(2:6:6*(N+1),1) = inf; %state y upper bound
obj.lbx(3:6:6*(N+1),1) = -inf; %state theta lower bound
obj.ubx(3:6:6*(N+1),1) = 6; %state theta upper bound
obj.lbx(4:6:6*(N+1),1) = -inf; %state x lower bound
obj.ubx(4:6:6*(N+1),1) = 10; %state x upper bound
obj.lbx(5:6:6*(N+1),1) = -inf; %state y lower bound
obj.ubx(5:6:6*(N+1),1) = 10; %state y upper bound
obj.lbx(6:6:6*(N+1),1) = -inf; %state theta lower bound
obj.ubx(6:6:6*(N+1),1) = 10; %state theta upper bound


obj.lbx(6*(N+1)+1:2:6*(N+1)+2*N,1) = v_f_min; %v lower bound
obj.ubx(6*(N+1)+1:2:6*(N+1)+2*N,1) = v_f_max; %v upper bound
obj.lbx(6*(N+1)+2:2:6*(N+1)+2*N,1) = v_b_min; %omega lower bound
obj.ubx(6*(N+1)+2:2:6*(N+1)+2*N,1) = v_b_max; %omega upper bound

        end

       %%  
       function u = stepImpl(obj,x,t)
disp(t)
N = 10;  
tic

xm = x(1:6);    % initial condition.
xr = x(7:12); % Reference posture.
u0 = zeros(N,2);        % two control inputs for each robot
X0 = repmat(xm,1,N+1)'; % initialization of the states decision variables
% initial value of the optimization variables
obj.x0  = [reshape(X0',6*(N+1),1);reshape(u0',2*N,1)];
obj.p   = [xm;xr]; % set the values of the parameters vector
solver = obj.casadi_solver;
sol = solver('x0', obj.x0, 'lbx', obj.lbx, 'ubx', obj.ubx,...
        'lbg', obj.lbg, 'ubg', obj.ubg,'p', obj.p);

u = reshape(full(sol.x(6*(N+1)+1:end))',2,N)';
u = u(1,:)';
toc
        end
%%
        function resetImpl(obj)
            % Initialize discrete-state properties.
        end
    end
end
