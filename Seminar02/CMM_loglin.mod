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
% This .mod file implements log-linearised 1st order conditions. The impulse 
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
c               ${c}$              (long_name='Log consumption')
ir              ${i}$              (long_name='Nominal interest rate')   
ir_ann          ${i^{ann}}$        (long_name='Nominal interest rate - annualised')  
n               ${n}$              (long_name='Log labour')
nu              ${\nu}$            (long_name='Log of monetary policy shock')
pie             ${\pi}$            (long_name='Inflation')   
pie_ann         ${\pi^{ann}}$      (long_name='Inflation') 
rr              ${r}$              (long_name='Real interest rate')  
rr_ann          ${r^{ann}}$        (long_name='Real interest rate - annualised')  
w               ${w}$              (long_name='Log real wage')
y               ${y}$              (long_name='Log output')
z               ${z}$              (long_name='Log of preference shifter')
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
gama          ${\gama}$            (long_name='Steady state growth')
phi           ${\phi}$             (long_name='Inverse of Frish elasticity') 
phi_pie       ${\phi_{\pi}}$       (long_name='Interest rate rule coefficient - inflation') 
rho           ${\rho}$             (long_name='Log inverse discount factor')   
rho_a         ${\rho_{a}}$         (long_name='AR persistence - technology') 
rho_nu        ${\rho_{\nu}}$       (long_name='AR persistence - monetary policy') 
rho_z         ${\rho_{z}}$         (long_name='AR persistence - preference shifter') 
siggma        ${\sigma}$           (long_name='Inverse of elasticity of substitution') 
siggma_a      ${\sigma_{a}}$       (long_name='Standard deviation - productivity shock') 
siggma_z      ${\sigma_{z}}$       (long_name='Standard deviation - preference shifter') 
siggma_nu     ${\sigma_{\nu}}$     (long_name='Standard deviation - monetary policy shock') 

pie_ss        ${\bar{\pi}}$        (long_name='Steady state - inflation') 
;


%% CALIBRATION (Chapter 3, page 52)
alppha    = 0.25;
betta     = 0.99;
gama      = 1;
phi       = 5;
phi_pie   = 1.5;
rho       = -log(betta);
rho_a     = 0.9;
rho_nu    = 0.5;
rho_z     = 0.5;
siggma    = 1;
siggma_a  = 1;
siggma_nu = 1;
siggma_z  = 1;

pie_ss    = 0;


%% MODEL
model(linear);

[name='Preference shifter (7)']
z = rho_z*z(-1) + e_z;

[name='Technology (11)']
a = rho_a*a(-1) + e_a;

[name='Consumption-Euler equation (26)']
c = c(+1) - 1/siggma*(ir - pie(+1) - rho) + 1/siggma*(1 - rho_z)*z;

[name='Labour supply (27)']
w = siggma*c + phi*n;

[name='Production function (28)']
y = a + (1 - alppha)*n;

[name='Labour demand (29)']
w = log(1 - alppha) + a - alppha*n;

[name='Real interest rate (37)']
rr = ir - pie(+1);

[name='Interest rate rule (43)']
ir = rho + pie_ss + phi_pie*(pie - pie_ss) + nu;

[name='Market clearing']
y = c;

[name='Monetary policy shock']
nu = rho_nu*nu(-1) + e_nu;

[name='Definition - annualised nominal interest rate']
ir_ann = 4*ir;

[name='Definition - annualised real interest rate']
rr_ann = 4*rr;

[name='Definition - annualised inflation']
pie_ann = 4*pie;

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
stoch_simul(order=1, irf=20) c y w n pie ir rr a;

