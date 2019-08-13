%% IDENTIFICATION
% JEM184 New Keynesian DSGE Modelling
% Institute of Economic Studies, Faculty of Social Sciences, Charles University
% Seminar 02b: A classical monetary model in Dynare
% Prepared by: Jan Zacek (jan.zacek@fsv.cuni.cz)
% Lastly modified: 13 August 2019


%% DESCRIPTION
% This file implements the basic Classical Monetary Model (CMM) presented 
% in Galí, J. (2015): Monetary Policy, Inflation, and the Business Cycle, 
% Princeton University Press, 2nd Edition, Chapter 2.
%
% This .mod file implements non-linear 1st order conditions. The impulse 
% response functions show log-linear deviations (percentage deviations)
% from steady state.


%% COPYRIGHT
% Copyright (C) 2019 Jan Zacek
% You can redistribute this software and/or modify it under the terms of
% the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.


%% VARIABLES
var 
a               ${a}$              (long_name='Log of technology')    
C               ${C}$              (long_name='Consumption')
N               ${N}$              (long_name='Labour')
PI              ${\Pi}$            (long_name='Gross inflation')   
Q               ${Q}$              (long_name='Price of bond')
IR              ${I}$              (long_name='Nominal interest rate')
RR              ${R}$              (long_name='Real interest rate')
Wr              ${W^{r}}$          (long_name='Real Wrage')
Y               ${Y}$              (long_name='Output')
z               ${z}$              (long_name='Log of preference shifter')
nu              ${\nu}$            (long_name='Log of monetary policy shock')
;


%% EXOGENOUS VARIABLES
varexo
e_a           ${\varepsilon^{a}}$  (long_name='Technology shock')
e_z           ${\varepsilon^{z}}$  (long_name='Preference shock')
e_nu          ${\varepsilon^{z}}$  (long_name='Monetary policy shock')
;


%% PARAMETERS
parameters
alppha        ${\alpha}$           (long_name='Share of labour')
betta         ${\beta}$            (long_name='Discount factor')
gama          ${\gama}$            (long_name='Steady state groWrth')
phi           ${\phi}$             (long_name='Inverse of Frish elasticity') 
phi_pi        ${\phi_{\pi}}$       (long_name='Interest rate rule coefficient - inflation') 
rho_a         ${\rho_{a}}$         (long_name='AR persistence - technology') 
rho_nu        ${\rho_{\nu}}$       (long_name='AR persistence - monetary policy') 
rho_z         ${\rho_{z}}$         (long_name='AR persistence - preference shifter') 
siggma        ${\sigma}$           (long_name='Inverse of elasticity of substitution') 
siggma_a      ${\sigma_{a}}$       (long_name='Standard deviation - productivity shock') 
siggma_z      ${\sigma_{z}}$       (long_name='Standard deviation - preference shifter') 
siggma_nu     ${\sigma_{\nu}}$     (long_name='Standard deviation - monetary policy shock') 

A_ss          ${\bar{A}}$
C_ss          ${\bar{C}}$
N_ss          ${\bar{N}}$
PI_ss         ${\bar{\Pi}}$
Q_ss          ${\bar{Q}}$
IR_ss
RR_ss
Wr_ss         ${\bar{W}^{r}}$
Y_ss          ${\bar{Y}}$
Z_ss          ${\bar{Z}}$
;


%% CALIBRATION (Chapter 3, page 52)
alppha    = 0.25;
betta     = 0.99;
gama      = 1;
phi       = 5;
phi_pi    = 1.5;
rho_a     = 0.9;
rho_nu    = 0.5;
rho_z     = 0.5;
siggma    = 1;
siggma_a  = 1;
siggma_z  = 1;
siggma_nu = 1;

A_ss      = 1;
N_ss      = (1 - alppha)^(1/(siggma*(1 - alppha) + alppha + phi));
Y_ss      = N_ss^(1-alppha);
C_ss      = Y_ss;
PI_ss     = 1;
Q_ss      = betta*gama^(-siggma)/PI_ss;
IR_ss     = 1/Q_ss;
RR_ss     = IR_ss/PI_ss;
Wr_ss     = (1-alppha)*N_ss^(-alppha);
Z_ss      = 1;


%% MODEL
model;

[name='Preference shifter (7)']
z = rho_z*z(-1) + e_z;

[name='Labour supply (8)']
exp(Wr) = exp(C)^(siggma)*exp(N)^(phi);

[name='Consumption-Euler equation (9)']
exp(Q) = betta*((exp(C(+1))/exp(C))^(-siggma)*(exp(z(+1))/(exp(z)))*1/exp(PI(+1)));

[name='Nominal interest rate (definition)']
exp(IR) = 1/exp(Q);

[name='Real interest rate (definition)']
exp(RR) = exp(IR)/exp(PI(+1));

[name='Production function (10)']
exp(Y) = exp(a)*exp(N)^(1-alppha);

[name='Technology (11)']
a = rho_a*a(-1) + e_a;

[name='Labour demand (13)']
exp(Wr) = (1-alppha)*exp(a)*exp(N)^(-alppha);

[name='Market clearing (15)']
exp(Y) = exp(C);

[name='Interest rate rule (43)']
exp(IR) = IR_ss*((exp(PI)/PI_ss)^(phi_pi))*exp(nu);

[name='Monetary policy shock']
nu = rho_nu*nu(-1) + e_nu;

end;


%% INITIAL VALUES
initval;
a  = log(A_ss);
C  = log(C_ss);
N  = log(N_ss);
PI = log(PI_ss);
Q  = log(Q_ss);
IR = log(IR_ss);
RR = log(RR_ss);
Wr = log(Wr_ss);
Y  = log(Y_ss);
z  = log(Z_ss);

end;


%% COVARIANCE MATRIX OF EXOGENOUS SHOCKS
shocks;
var e_a  = siggma_a^2;
var e_z  = siggma_z^2;
var e_nu = siggma_nu^2;
end;


%% STEADY STATE COMPUTATION AND STABILITY CHECK
steady;
check;
resid(1);


%% LATEX OUTPUT COMMANDS
write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;
collect_latex_files;


%% STOCHASTIC SIMULATION 
stoch_simul(order=1, irf=20) C Y Wr N PI IR RR a;

