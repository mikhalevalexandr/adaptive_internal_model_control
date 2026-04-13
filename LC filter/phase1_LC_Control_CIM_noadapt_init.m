
% clear; clc; close all;

% Parameters (Table 3.1)
L1 = 150e-6;     % H inverter side inductance
R1 = 0.045;%0.1;      % ohm
Rd = 1;        % ohm damping resistance
Cf = 22e-6;      % F

% % State matrix A (3x3)
A = [ -(R1+Rd)/L1,       -1/L1;
        1/Cf,            0    ];

% B1 for w=[v_ref; v_g] (3x2)
B1 = [  1/L1;
        0];

% B2 for u (3x1)
B2 = [ 1/L1;
        0 ];

% Output y = e (1x3), with C1, D1, D2 as in the screenshot
C = [1, 0;
     0, 1];     % C1
D1 = [0;
      0];        % maps w to y
D2 = [0;
      0];             % maps u to y

% Build ss with combined input [w; u] = [v_ref; v_g; u]
B = [B1, B2];       % 3x3
D = [D1, D2];       % 1x3
P = ss(A, B, C, D);

s = tf('s');

%% IMC
Pu = ss(A, B2, C, D2);
A=A;
B=B2;
C=C;
D=D2;

Putf = tf(Pu);
p_des = [-50000, -10000];

K = -place(A,B2,p_des);
Acl = A + B2*K;
Bcl = B2;          % for exogenous channels we’ll add in the 2-DOF section
Ccl = C;
Dcl = D2;
eigAcl = eig(Acl)


E  = 1;         % (p_e = 1), measured y has p = 2
Ecomp = 0;      % E? (complement) – any matrix s.t. [E;Ecomp] is orthonormal-ish row-wise

p_e=size(E,1);
p=size(C,1);
n=size(A,1);
m=size(B,2);
% Internal model M(s) 1 freq
freq=1;
w=2*pi*freq;

am=2*pi;
kef=8;
% phi=(s+am*kef)^4;
phi=(s+am*kef)^2;
w1=w;
w2=3*w;

Mtf = phi / ( (s^2 + w2^2) );
[Am,Bm,Cm,Dm] = ssdata(ss(Mtf));
M  = ss(Am,Bm,Cm,Dm);

nu = size(Bm,2);    % number of parameters
n_m=size(Am,1);

Bsharp = (Pu.B' * Pu.B)\Pu.B';     

A_c=Am;
C_c=Cm;
B_c=[Bm Bm*(K+Bsharp*A)-Am*Bm*Bsharp];
D_c=[eye(nu) K-Cm*Bm*Bsharp];

B_c_st=Bm*(K+Bsharp*A)-Am*Bm*Bsharp;
D_c_st=K-Cm*Bm*Bsharp;
SFControl=ss(A_c, B_c, C_c, D_c);
SFControl_st=ss(A_c, B_c_st, C_c, D_c_st);
Pi=(1-inv(M))*(Bsharp*s-Bsharp*A);
Pimin=minreal(Pi,1e-6);
Piminss=dss2ss(Pimin);
Piminssmin =minreal(Piminss,1e-6);
Pitf=zpk(Piminssmin);
v2u=Mtf/(1-Mtf*(K-Pi)*Putf);
dcg=dcgain(v2u);
dcgain(1/(1-K*Putf))
v2umin=minreal(v2u,1e-6);


freq0var=freq*3;