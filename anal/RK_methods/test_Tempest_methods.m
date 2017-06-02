% Simple script to instantiate and test the order conditions for
% the RK methods in Tempest+ARKode
%
% Daniel R. Reynolds
% SMU Mathematics
% August 2016

clear

% KGU(6,3)
s = 6;
q = 3;
A = zeros(s,s);
A(2,1) = 0.2;
A(3,2) = 0.2;
A(4,3) = 1/3;
A(5,4) = 2/3;
A(6,1) = 0.25;  A(6,5) = 0.75;
c = [0; 0.2; 0.2; 1/3; 2/3; 1];
b = [0.25, 0, 0, 0, 0.75, 0];
fprintf('\n\nChecking KGU(6,3) via check_erk:\n');
[q,lq] = check_erk(c,A,b,1e-11,1,true,[-4,1,-5,5],'KGU(6,3)','kgu63');
!epstopdf kgu63_stab_region.eps
!rm kgu63_stab_region.eps
[qse,qsi,qsa] = check_ark_stage_order(c,c,A,A,1e-11);
fprintf('  stage order = %i\n', qse);
fprintf('  max stable step on imaginary axis = %g\n', ...
        max_stable_step(A,b,sqrt(-1)));


% SSPRK(5,4)
s = 5;
q = 4;
alpha1 = 0.555629506348765;
alpha2 = 0.379898148511597;
alpha3 = 0.821920045606868;
beta1 = 0.517231671970585;
beta2 = 0.096059710526147;
beta3 = 0.386708617503259;
A = zeros(s,s);
A(2,1) = 0.391752226571890; 
A(3,1:2) = [alpha1*A(2,1), 0.368410593050371];
A(4,1:3) = [alpha2*A(3,1:2), 0.251891774271694];
A(5,1:4) = [alpha3*A(4,1:3), 0.544974750228521];
c = zeros(s,1);      
c(2) = A(2,1);
c(3) = sum(A(3,:));
c(4) = sum(A(4,:));
c(5) = sum(A(5,:));
b = zeros(1,s);
b(1) = beta1*A(3,1) + beta2*A(4,1) + beta3*A(5,1);
b(2) = beta1*A(3,2) + beta2*A(4,2) + beta3*A(5,2);
b(3) = beta2*A(4,3) + beta3*A(5,3);
b(4) = beta3*A(5,4) + 0.063692468666290;
b(5) = 0.226007483236906;
fprintf('\n\nChecking SSPRK(5,4) via check_erk:\n');
[q,lq] = check_erk(c,A,b,1e-11,1,true,[-6,1,-5,5],'SSPRK(5,4)','ssprk54');
!epstopdf ssprk54_stab_region.eps
!rm ssprk54_stab_region.eps
[qse,qsi,qsa] = check_ark_stage_order(c,c,A,A,1e-11);
fprintf('  stage order = %i\n', qse);
fprintf('  max stable step on imaginary axis = %g\n', ...
        max_stable_step(A,b,sqrt(-1)));


% ARS233
s = 3;
q = 3;
gamma = (3 + sqrt(3)) / 6;
ci = [0; gamma; 1.0 - gamma];
Ai = zeros(s,s);
Ai(2,2) = gamma;
Ai(3,2:3) = [1 - 2 * gamma, gamma];
bi = [0, 0.5, 0.5];     
ce = ci;
Ae = zeros(s,s);
Ae(2,1) = gamma;
Ae(3,1:2) = [gamma - 1; 2*(1 - gamma)];
be = bi;
fprintf('\n\nChecking ARS233 via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-5,2,-4,4],'ARS233','ars233');
!epstopdf ars233_stab_regions.eps
!rm ars233_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));


% ARS232
s = 3;
q = 2;
gamma = 1 - 1/sqrt(2);
delta = -2 * sqrt(2) / 3;
ci = [0.0; gamma; 1.0];
Ai = zeros(s,s);      
Ai(2,2) = gamma;
Ai(3,2:3) = [1 - gamma, gamma];
bi = [0, 1 - gamma, gamma];
ce = ci;
Ae = zeros(s,s);
Ae(2,1) = gamma;
Ae(3,1:2) = [delta, 1 - delta];
be = bi;      
fprintf('\n\nChecking ARS232 via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-4,2,-4,4],'ARS232','ars232');
!epstopdf ars232_stab_regions.eps
!rm ars232_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));


% ARS222
s = 3;
q = 2;
gamma = 1 - 1/sqrt(2);
delta = 1 - 1/(2*gamma);
ci = [0; gamma; 1];
Ai = zeros(s,s);      
Ai(2,2) = gamma;
Ai(3,2:3) = [1 - gamma, gamma];
bi = [0, 1 - gamma, gamma];
ce = ci;
Ae = zeros(s,s);
Ae(2,1) = gamma;
Ae(3,1:2) = [delta, 1 - delta];
be = [delta, 1 - delta, 0];
fprintf('\n\nChecking ARS222 via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-4,2,-3,3],'ARS222','ars222');
!epstopdf ars222_stab_regions.eps
!rm ars222_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));


% ARS343
s = 4;
q = 3;
gamma  = 0.4358665215084590;
gamma2 = gamma * gamma;
b1 = -1.5 * gamma2 + 4.0 * gamma - 0.25;
b2 =  1.5 * gamma2 - 5.0 * gamma + 1.25;
a42 = 0.5529291480359398;
a43 = 0.5529291480359398;
a31 = (1.0 - 4.5 * gamma + 1.5 * gamma2) * a42 + (2.75 - 10.5 * gamma + 3.75 * gamma2) * a43 - 3.5 + 13 * gamma - 4.5 * gamma2;
a32 = (-1.0 + 4.5 * gamma - 1.5 * gamma2) * a42	+ (-2.75 + 10.5 * gamma - 3.75 * gamma2) * a43 + 4.0 - 12.5 * gamma + 4.5 * gamma2;
a41 = 1.0 - a42 - a43;
ci = [0; gamma; (1+gamma)/2; 1];
Ai = [0, 0, 0, 0; 0, gamma, 0, 0; 0, (1-gamma)/2, gamma, 0; 0, b1, b2, gamma]; 
bi = Ai(4,:);
ce = ci;
Ae = [0, 0, 0, 0; gamma, 0, 0, 0; a31, a32, 0, 0; a41, a42, a43, 0];
be = bi;
fprintf('\n\nChecking ARS343 via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-5,2,-5,5],'ARS343','ars343');
!epstopdf ars343_stab_regions.eps
!rm ars343_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% ARS443
s = 5;
q = 3;
ci = [0; 0.5; 2/3; 0.5; 1];
Ai = [0, 0, 0, 0, 0; 
      0, 0.5, 0, 0, 0; 
      0, 1/6, 0.5, 0, 0;
      0, -0.5, 0.5, 0.5, 0;
      0, 1.5, -1.5, 0.5, 0.5];
bi = Ai(5,:);      
ce = ci;
Ae = [0, 0, 0, 0, 0; 
      0.5, 0, 0, 0, 0;
      11/18, 1/18, 0, 0, 0;
      5/6, -5/6, 0.5, 0, 0;
      0.25, 7/4, 0.75, -7/4, 0];
be = Ae(5,:);
fprintf('\n\nChecking ARS443 via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-4,2,-4,4],'ARS443','ars443');
!epstopdf ars443_stab_regions.eps
!rm ars443_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% ARK232
s = 3;
q = 2;
gamma = 1 - 1/sqrt(2);
alpha = 1/6 * (3 + 2 * sqrt(2));
delta = 1 / (2 * sqrt(2));
twogamma = 2 * gamma;
ci = [0; twogamma; 1];
Ai = [0, 0, 0; gamma, gamma, 0; delta, delta, gamma];
bi = Ai(3,:);
ce = ci;
Ae = [0, 0, 0; twogamma, 0, 0; 1 - alpha, alpha, 0];
be = bi;      
fprintf('\n\nChecking ARK232 via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-5,2,-5,5],'ARK232','ark232');
!epstopdf ark232_stab_regions.eps
!rm ark232_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% SSP2(2,2,2)
s = 2;
q = 2;
gamma = 1 - 1 / sqrt(2);
ci = [gamma; 1.0 - gamma];
Ai = [gamma, 0; 1 - 2 * gamma, gamma];
bi = 0.5*[1 1];
ce = [0; 1];
Ae = [0, 0; 1, 0];
be = bi;
fprintf('\n\nChecking SSP2(2,2,2) via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-4,2,-5,5],'SSP2(2,2,2)','ssp2_222');
!epstopdf ssp2_222_stab_regions.eps
!rm ssp2_222_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% SSP2(3,3,2)a
s = 3;
q = 2;
ci = [0.25; 0.25; 1.0];
Ai = [0.25, 0, 0; 0, 0.25, 0; 1/3, 1/3, 1/3];
bi = Ai(3,:);
ce = [0; 0.5; 1.0];
Ae = [0, 0, 0; 0.5, 0, 0; 0.5, 0.5, 0];
be = bi;
fprintf('\n\nChecking SSP2(3,3,2)a via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-6,1,-5,5],'SSP2(3,3,2)a','ssp2_332a');
!epstopdf ssp2_332a_stab_regions.eps
!rm ssp2_332a_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));




% SSP3(3,3,2)
s = 3;
q = 2;
gamma = 1 - 1 / sqrt(2); 
ci = [gamma; 1 - gamma; 0.5];
Ai = [gamma, 0, 0; 1 - 2 * gamma, gamma, 0; 0.5 - gamma, 0, gamma];
bi = [1 / 6, 1 / 6, 2 / 3];
ce = [0; 1; 0.5];
Ae = [0, 0, 0; 1, 0, 0; 0.25, 0.25, 0];
be = bi;
fprintf('\n\nChecking SSP3(3,3,2) via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-4,2,-4,4],'SSP3(3,3,2)','ssp3_332');
!epstopdf ssp3_332_stab_regions.eps
!rm ssp3_332_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% SSP3(4,3,3)
s = 4;
q = 3;
alpha = 0.24169426078821; 
beta  = 0.06042356519705;
eta   = 0.12915286960590;
delta = 0.5 - beta - eta - alpha;    
ci = [alpha; 0.0; 1; 0.5];
Ai = zeros(s,s);      
Ai = [alpha, 0, 0, 0; -alpha, alpha, 0, 0; 0, 1 - alpha, alpha, 0; beta, eta, delta, alpha]; 
bi = [0, 1/6, 1/6, 2/3];
ce = [0; 0; 1; 0.5];
Ae = [0, 0, 0, 0; 0, 0, 0, 0; 0, 1, 0, 0; 0 0.25, 0.25, 0];
be = bi;
fprintf('\n\nChecking SSP3(4,3,3) via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-4,2,-5,5],'SSP3(4,3,3)','ssp3_433');
!epstopdf ssp3_433_stab_regions.eps
!rm ssp3_433_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% SSP2(3,3,2)b
s = 3;
q = 2;
ci = [0.2; 0.3; 1.0];
Ai = zeros(s,s);      
Ai(1,1) = 0.2;
Ai(2,1:2) = [0.1, 0.2];
Ai(3,:) = [1/3, 1/3, 1/3];
bi = Ai(3,:);
ce = [0; 0.5; 1];
Ae = zeros(s,s);
Ae(2,1) = 0.5;
Ae(3,1:2) = [0.5, 0.5];
be = bi;
fprintf('\n\nChecking SSP2(3,3,2)b via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-6,1,-4,4],'SSP2(3,3,2)b','ssp2_332b');
!epstopdf ssp2_332b_stab_regions.eps
!rm ssp2_332b_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% SSP3(3,3,3)
s = 3;
q = 3;
ci = [0; 1; 0.5];
Ai = [0, 0, 0; 14/15, 1/15, 0; 7/30, 0.2, 1/15]; 
bi = [1/6, 1/6, 2/3];
ce = ci;
Ae = [0, 0, 0; 1, 0, 0; 0.25, 0.25, 0];      
be = bi;
fprintf('\n\nChecking SSP3(3,3,3) via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-6,1,-4,4],'SSP3(3,3,3)','ssp3_333');
!epstopdf ssp3_333_stab_regions.eps
!rm ssp3_333_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% SSP2(3,3,2)-LSPUM
s = 3;
q = 2;
ci = [2/11; 289/462; 751/924];
Ai = [2/11, 0, 0; 205/462, 2/11, 0; 2033/4620, 21/110, 2/11]; 
bi = [24/55, 0.2, 4/11];
ce = [0; 5/6; 11/12];
Ae = [0, 0, 0; 5/6, 0, 0; 11/24, 11/24, 0];
be = bi;
fprintf('\n\nChecking SSP2(3,3,2)-LSPUM via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-5,2,-5,5],'SSP2(3,3,2)-LSPUM','ssp2_332_lpsum');
!epstopdf ssp2_332_lpsum_stab_regions.eps
!rm ssp2_332_lpsum_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% SSP2(3,3,2)-LPUM
s = 3;
q = 2;
ci = [2/11; 69/154; 67/77];
Ai = [2/11, 0, 0; 41/154, 2/11, 0; 289/847, 42/121, 2/11]; 
bi = 1/3*[1, 1, 1];
ce = [0; 0.5; 1.0];
Ae = [0, 0, 0; 0.5, 0, 0; 0.5, 0.5, 0];
be = bi;
fprintf('\n\nChecking SSP2(3,3,2)-LPUM via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-6,1,-5,5],'SSP2(3,3,2)-LPUM','ssp2_332_lpum');
!epstopdf ssp2_332_lpum_stab_regions.eps
!rm ssp2_332_lpum_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% SSP2(3,3,2)-LPM1
s = 3;
q = 2;
ci = [2/11; 4523/9317; 15517/18634];
Ai = [2/11, 0, 0; 2829/9317, 2/11, 0; 148529/428582, 7/23, 2/11]; 
bi = 1/3*[1, 1, 1];
ce = [0; 0.5; 1];
Ae = [0, 0, 0; 0.5, 0, 0; 0.5, 0.5, 0];
be = bi;
fprintf('\n\nChecking SSP2(3,3,2)-LPM1 via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-6,1,-5,5],'SSP2(3,3,2)-LPM1','ssp2_332_lpm1');
!epstopdf ssp2_332_lpm1_stab_regions.eps
!rm ssp2_332_lpm1_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% SSP2(3,3,2)-LPM2
s = 3;
q = 2;
ci = [2/11; 5003/13310; 6271/6655];
Ai = [2/11, 0, 0; 2583/13310, 2/11, 0; 39731/139755, 10/21, 2/11]; 
bi = 1/3*[1, 1, 1];
ce = [0; 0.5; 1];
Ae = [0, 0, 0; 0.5, 0, 0; 0.5, 0.5, 0];      
be = bi;
fprintf('\n\nChecking SSP2(3,3,2)-LPM2 via check_ark2:\n');
[qE,qI,q] = check_ark2(ce,ci,Ae,Ai,be,bi,1e-11,1,true,[-6,1,-5,5],'SSP2(3,3,2)-LPM2','ssp2_332_lpm2');
!epstopdf ssp2_332_lpm2_stab_regions.eps
!rm ssp2_332_lpm2_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% ARK3(2)4L[2]SA
B = butcher('ARK3(2)4L[2]SA-ERK');
[m,n] = size(B);
s = n-1;
ce = B(1:s,1);
Ae = B(1:s,2:s+1);
be = B(s+1,2:s+1);
de = B(s+2,2:s+1);
B = butcher('ARK3(2)4L[2]SA-ESDIRK');
ci = B(1:s,1);
Ai = B(1:s,2:s+1);
bi = B(s+1,2:s+1);
di = B(s+2,2:s+1);
fprintf('\n\nChecking ARK3(2)4L[2]SA via check_ark_embedded2:\n');
[qE,qI,q,pE,pI,p] = check_ark_embedded2(ce,ci,Ae,Ai,be,bi,de,di,1e-11,1,true,[-5,2,-5,5],'ARK3(2)4L[2]SA','ark324l2sa');
!epstopdf ark324l2sa_stab_regions.eps
!rm ark324l2sa_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% ARK4(3)6L[2]SA
B = butcher('ARK4(3)6L[2]SA-ERK');
[m,n] = size(B);
s = n-1;
ce = B(1:s,1);
Ae = B(1:s,2:s+1);
be = B(s+1,2:s+1);
de = B(s+2,2:s+1);
B = butcher('ARK4(3)6L[2]SA-ESDIRK');
ci = B(1:s,1);
Ai = B(1:s,2:s+1);
bi = B(s+1,2:s+1);
di = B(s+2,2:s+1);
fprintf('\n\nChecking ARK4(3)6L[2]SA via check_ark_embedded2:\n');
[qE,qI,q,pE,pI,p] = check_ark_embedded2(ce,ci,Ae,Ai,be,bi,de,di,1e-11,1,true,[-6,2,-6,6],'ARK4(3)6L[2]SA','ark436l2sa');
!epstopdf ark436l2sa_stab_regions.eps
!rm ark436l2sa_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));



% ARK5(4)8L[2]SA
B = butcher('ARK5(4)8L[2]SA-ERK');
[m,n] = size(B);
s = n-1;
ce = B(1:s,1);
Ae = B(1:s,2:s+1);
be = B(s+1,2:s+1);
de = B(s+2,2:s+1);
B = butcher('ARK5(4)8L[2]SA-ESDIRK');
ci = B(1:s,1);
Ai = B(1:s,2:s+1);
bi = B(s+1,2:s+1);
di = B(s+2,2:s+1);
fprintf('\n\nChecking ARK5(4)8L[2]SA via check_ark_embedded2:\n');
[qE,qI,q,pE,pI,p] = check_ark_embedded2(ce,ci,Ae,Ai,be,bi,de,di,1e-11,1,true,[-8,0.5,-7,7],'ARK5(4)8L[2]SA','ark548l2sa');
!epstopdf ark548l2sa_stab_regions.eps
!rm ark548l2sa_stab_regions.eps
[qse,qsi,qsa] = check_ark_stage_order(ce,ci,Ae,Ai,1e-11);
fprintf('  stage order: ERK = %i, DIRK = %i, ARK = %i\n', qse, qsi, qsa);
fprintf('  max stable explicit step on imaginary axis = %g\n', ...
        max_stable_step(Ae,be,sqrt(-1)));
