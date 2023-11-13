% Clear contents
clear; clc

% Create symbolic variables
syms t1(t) t2(t) t4(t) t3(t) %H1 H2 L1 L2;

% Define known variables
H1 = 89.45;
H2 = 100;
L1 = 35;
L2 = 100;

% Solve for the Jacobian of joint 1
sw1  = [0; 0; 1];
a1   = [0; 0; 0];
sv1  = cross(a1, sw1);
S1   = [sw1; sv1];
sw1_brack = [      0, -sw1(3),  sw1(2);
              sw1(3),       0, -sw1(1);
             -sw1(2),  sw1(1),       0];
J1 = S1;

% Solve for the Jacobian of joint 2
sw2       = [0; 1; 0];
a2        = [0; 0; H1];
sv2       = cross(a2, sw2);
S2        = [sw2; sv2];
sw2_brack = [      0, -sw2(3),  sw2(2);
              sw2(3),       0, -sw2(1);
             -sw2(2),  sw2(1),       0];

R1       = eye(3) + sin(t1)*sw1_brack + (1-cos(t1))*(sw1_brack*sw1_brack);
p1       = (eye(3)*t1 + (1-cos(t1))*sw1_brack + (t1-sin(t1))*(sw1_brack*sw1_brack))*sv1;
p1       = p1(t);
p1_brack = [     0, -p1(3),  p1(2);
             p1(3),      0, -p1(1);
            -p1(2),  p1(1),     0];
adj1     = horzcat(R1, zeros(3));
adj1     = vertcat(adj1, [p1_brack*R1, R1]);
J2       = adj1 * S2;

% Solve for the Jacobian of joint 3
sw3       = [0; 1; 0];
a3        = [L1; 0; H1+H2];
sv3       = cross(a3, sw3);
S3        = [sw3; sv3];
sw3_brack = [      0, -sw3(3),  sw3(2);
              sw3(3),       0, -sw3(1);
             -sw3(2),  sw3(1),       0];

R2       = eye(3) + sin(t2)*sw2_brack + (1-cos(t2))*(sw2_brack*sw2_brack);
p2       = (eye(3)*t2 + (1-cos(t2))*sw2_brack + (t2-sin(t2))*(sw2_brack*sw2_brack))*sv2;
p2       = p2(t);
p2_brack = [     0, -p2(3),  p2(2);
             p2(3),      0, -p2(1);
            -p2(2),  p2(1),     0];
adj2     = horzcat(R2, zeros(3));
adj2     = vertcat(adj2, [p2_brack*R2, R2]);
J3       = adj2 * S3;

% Solve for the Jacobian of joint 4
sw4       = [0; 1; 0];
a4        = [L1+L2; 0; H1+H2];
sv4       = cross(a4, sw4);
S4        = [sw4; sv4];
sw4_brack = [      0, -sw4(3),  sw4(2);
              sw4(3),       0, -sw4(1);
             -sw4(2),  sw4(1),       0];

R3       = eye(3) + sin(t3)*sw3_brack + (1-cos(t3))*(sw3_brack*sw3_brack);
p3       = (eye(3)*t3 + (1-cos(t3))*sw3_brack + (t3-sin(t3))*(sw3_brack*sw3_brack))*sv3;
p3       = p3(t);
p3_brack = [     0, -p3(3),  p3(2);
             p3(3),      0, -p3(1);
            -p3(2),  p3(1),     0];
adj3     = horzcat(R3, zeros(3));
adj3     = vertcat(adj3, [p3_brack*R3, R3]);
J4       = adj3 * S4;

% Concatenate Jacobians of joints 1 and 2
J = [J1, J2, J3, J4]
J = vpa(J)

% Solve for Vs
q = [t1; t2; t3; t4];
q_dot = simplify(diff(q, t));

disp("Second solution of Vs:")
Vs2 = J * q_dot