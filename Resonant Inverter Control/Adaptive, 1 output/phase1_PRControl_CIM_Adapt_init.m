clear; clc; close all;

% Parameters
L1 = 150e-6;     % H inverter side inductance
R1 = 0.045;      % ohm
Rd = 1;          % ohm damping resistance
Cf = 22e-6;      % F
L2 = 450e-6;     % H grid side inductance
R2 = 0.135;      % ohm
Lg=0;            % H line inductance of the grid
Rg=0;
L2g=L2+Lg;
R2g=R2+Rg;

A = [ -(R1+Rd)/L1,     Rd/L1,       -1/L1;
        Rd/L2g,     -(R2g+Rd)/L2g,      1/L2g;
        1/Cf,        -1/Cf,         0    ];


B1 = [  1/L1;
       -1/L2g;
        0];

B2 = [ 1/L1;
        0;
        0 ];

C = [0, 1, 0]; 
D1 = 0;        
D2 = 0;        

B = [B1, B2];       
D = [D1, D2];       
P = ss(A, B, C, D);



s = tf('s');
Gfffilter=33*(0.05*s+1)/((s+300)*(0.002*s+1));%feed-forward filter

%% CIM
Pu = ss(A, B2, C, D2);
Pu_tf=tf(Pu);
A=A;
B=B2;
C=C;
D=D2;

E  = 1;        
Ecomp = 0;      

p_e=size(E,1);
p=size(C,1);
n=size(A,1);
m=size(B,2);

freq=49.5;
w=2*pi*freq;

am=2*pi*80;


w1=w;
w2=3*w;
w3=5*w;
freq0=50;
w0=2*pi*freq0;
w01=w0;
w02=3*w0;
w03=5*w0;


%% 1 harm model generation

% 1 harm

phi=(s+am)^2;
Theta=w1^2;
Theta0=w0^2;

M0_tf=s^2/phi;
M1_tf=1/phi;




sys0 = ss(M0_tf);
sys1 = ss(M1_tf);

% dual systems
sys0d = ss(sys0.A', sys0.C', sys0.B', sys0.D');
sys1d = ss(sys1.A', sys1.C', sys1.B', sys1.D');

% controllable companion on dual
sys0d_ccf = canon(sys0d,'companion');
sys1d_ccf = canon(sys1d,'companion');

% transpose back -> observable companion
sys0_ocf = ss(sys0d_ccf.A', sys0d_ccf.C', sys0d_ccf.B', sys0d_ccf.D');
sys1_ocf = ss(sys1d_ccf.A', sys1d_ccf.C', sys1d_ccf.B', sys1d_ccf.D');

Am0_can = sys0_ocf.A;
Bm0_can = sys0_ocf.B;
Cm0_can = sys0_ocf.C;
Dm0_can = sys0_ocf.D;
M0_can=ss(Am0_can,Bm0_can,Cm0_can,Dm0_can);
M0_can_tf=tf(M0_can);
Am1_can = sys1_ocf.A;
Bm1_can = sys1_ocf.B;
Cm1_can = sys1_ocf.C;
Dm1_can = sys1_ocf.D;
M1_can=ss(Am0_can,Bm1_can,Cm0_can,Dm1_can);
M1_can_tf=tf(M1_can);








n_m0=size(Am0_can,1);
m_m1=size(Bm1_can,2);

Bsharp = (B' * B)\B';     % m x n  
N = null(B');             % columns span orthogonal complement
Bperp = N.';              % (n-m) x n
% Check nonsingularity of [Bperp; B#]:
if abs(det([Bperp; Bsharp])) < 1e-10
    error('Singularity in [B?; B#]. Choose a different basis for B?.');
end
% Solve generalized Sylvester for X

Asylv=[Bperp;zeros(p_e,n)];
Bsylv=Am0_can;

Csylv=-[Bperp*A;-E*C];
Dsylv=eye(size(Bsylv));
Esylv=[zeros(n-m,n_m0);-Cm0_can];
X = generalized_sylvester(Asylv, Bsylv, Csylv, Dsylv, Esylv);


A_np = A + X*Bm0_can*E*C;
B_np = B;
C_np = C;
D_np = D;
Pubar = ss(A_np, B_np, C_np, D_np);
Pubartf=tf(Pubar);
Puzt = tzero(Pu);     % transmission zeros
Pupt = pole(Pu);      % poles of the realization
C0_hat = Bsharp * X * (Am0_can) - Bsharp * A * X;


Cpai_hat = Bsharp * X * (Am0_can) - Bsharp * A * X;



%Compensators
Pai1_hat=ss(Am0_can, Bm0_can, Cpai_hat, zeros(m,p_e)); 

Cpai2_hat = Ecomp*C*X+Ecomp*D*Cpai_hat;

Pai2_hat=ss(Am0_can, Bm0_can, Cpai2_hat, zeros(p_e));

Rubar=-5;
R_ss = ss([],[],[], Rubar);
R_tf=tf(R_ss);
RI_tf=[eye(p_e) R_tf];

%% Psi0 in SS

Edss = blkdiag( eye(n_m0+n_m0+n),zeros(p_e));

Adss = [ Am0_can-Bm0_can*Cm0_can,  zeros(n_m0),              zeros(n_m0,n), zeros(n_m0,m);
              zeros(n_m0),             Am0_can-Bm0_can*Cm0_can,  zeros(n_m0,n), zeros(n_m0,m);
              zeros(n,n_m0)            zeros(n,n_m0),              A_np,            B;
              zeros(p_e,n_m0)          -Cm0_can,                   -E*C,           zeros(p_e,m)];


Bdss = [ Bm0_can*E zeros(n_m0,m_m1);
         zeros(n_m0,p) -Bm1_can;
         zeros(n,p) zeros(n,m_m1);
         zeros(p_e,p) zeros(p_e,m_m1)];

Cdss = [ -Cpai_hat, zeros(p_e,n_m0), zeros(p_e,n), eye(p_e,m);
        Ecomp'*Cpai2_hat-E'*Cm0_can, zeros(p,n_m0),     -C,             -D;
        -Cm0_can, Cm0_can, zeros(p_e,n), zeros(p_e,m)];

Ddss = [zeros(p_e,p) zeros(p_e,m_m1);
        eye(p) zeros(p,m_m1);
        E      zeros(p_e,m_m1)];
Psi0dss_dss=dss(Adss,Bdss,Cdss,Ddss,Edss);
Psi0dssss=dss2ss(Psi0dss_dss);
Psi0dssss_min=minreal(Psi0dssss, 1e-8);
Psi0dssss_tf=tf(Psi0dssss);


% Adaptation params
Q_adapt=zeros(m_m1,1); alpha_adapt=100; k0_adapt=1000; c_adapt=0.1;
PiInit        = k0_adapt*eye(m_m1);