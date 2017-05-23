function [q,lq] = check_erk(c,A,b,tol,reportL,doPlot,box,mname,fname)
% Usage: [q,lq] = check_erk(c,A,b,tol,reportL,doPlot,box,mname,fname)
% 
% Checks the explicit RK table given by
%
%     c | A
%    -------
%       | b
%
% to determine the analytical order of accuracy (up to 6th).  We 
% separately check for "linear order of accuracy" on the linear
% autonomous ODE 
%    y' = A*y
% and for actual order of accuracy, using analytical order
% conditions from Sandu & Gunther, SINUM 53, 2015.
%
% Addtionally, this routine additionally checks for added stability
% conditions for ERK methods [Kinnmark & Gray, 1984]
%
% Inputs:
%   c, A, b -- ERK table
%   tol -- tolerance for checking order conditions (use negative
%          or zero for default value)
%   reportL -- integer flag denoting print level:
%              >1 -> all results
%              =1 -> final results
%            else -> no results
%   doPlot -- boolean flag denoting whether to create/save plot
%   box -- [xl,xr,yl,yr] sub-region of the complex plane to
%          use in plot (ignored if doPlot is false)
%   mname -- string containing method name to insert into plot title
%   fname -- string containing method name to insert into plot filenames
%
% Outputs:
%   q -- order of accuracy for ERK method
%   lq -- linear order of accuracy for ERK method
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2016, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------


% set tolerance on 'equality'
if (tol <= 0)
   tol = 1e-8;
end

% initialize failure flags, order of accuracy
Ofail = false;
Lfail = false;
q = -1;
ql = -1;

% get number of stages
s = length(b);

% convert b to column vector for these tests
b = reshape(b,s,1);

% check row sum condition
for i=1:s
   tst = c(i) - sum(A(i,:));
   if (abs(tst) > tol)
      Ofail = true;
      Lfail = true;
      if (reportL>1)
         fprintf('    Method fails row sum condition, i = %i, tst = %g\n', i, tst);
      end
   end
end
if (~Ofail)
  q = 0;
  lq = 0;
end

% check for first order
if (~Ofail || ~Lfail)
   tst = sum(b)-1;
   if (abs(tst) > tol)
      Ofail = true;
      Lfail = true;
      if (reportL>1)
         fprintf('    Method fails 1st order condition (tst = %g)\n', tst);
      end
   end
   if (reportL>1)
      if (~Ofail), fprintf('  Method passes order 1 conditions\n'); end
   end
   if (~Ofail)
      q = 1;
      lq = 1;
   end
end

% check for second order
if (~Ofail || ~Lfail)
   tst = b'*c - 0.5;
   if (abs(tst) > tol)
      Ofail = true;
      Lfail = true;
      if (reportL>1)
         fprintf('    Method fails 2nd order condition (tst = %g)\n', tst);
      end
   end
   if (reportL>1)
      if (~Ofail),  fprintf('  Method passes order 2 condition\n'); end
   end
   if (~Ofail)
      q = 2;
      lq = 2;
   end
end

% check for third order
if (~Ofail || ~Lfail) 
   
   tst = b'*(c.*c) - (1/3);
   if (abs(tst) > tol)
      Ofail = true;
      if (reportL>1)
         fprintf('    Method fails 3nd order condition (tst = %g)\n', tst);
      end
   end

   tst = b'*(A*c) - (1/6);
   if (abs(tst) > tol)
      Ofail = true;
      Lfail = true;
      if (reportL>1)
         fprintf('    Method fails linear 3nd order condition (tst = %g)\n', tst);
      end
   end

   if (reportL>1)
      if (~Ofail),  fprintf('  Method passes order 3 conditions\n'); end
      if (~Lfail && Ofail),  fprintf('  Method passes linear order 3 condition\n'); end
   end
   if (~Ofail)
      q = 3;
   end
   if (~Lfail)
      lq = 3;
   end
end

% check for third order
if (~Ofail || ~Lfail) 

   if (~Ofail)
      tst = b'*(c.*c.*c) - (1/4);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 4th order condition A (tst = %g)\n', tst);
         end
      end

      tst = (b.*c)'*(A*c) - (1/8);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 4th order condition B (tst = %g)\n', tst);
         end
      end

      tst = b'*A*(c.*c) - (1/12);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 4th order condition C (tst = %g)\n', tst);
         end
      end
   end

   tst = b'*A*A*c - (1/24);
   if (abs(tst) > tol)
      Ofail = true;
      Lfail = true;
      if (reportL>1)
         fprintf('    Method fails linear 4th order condition (tst = %g)\n', tst);
      end
   end

   if (reportL>1)
      if (~Ofail),  fprintf('  Method passes order 4 conditions\n'); end
      if (~Lfail && Ofail),  fprintf('  Method passes linear order 4 condition\n'); end
   end
   if (~Ofail)
      q = 4;
   end
   if (~Lfail)
      lq = 4;
   end
end

% check for fifth order
if (~Ofail || ~Lfail) 

   if (~Ofail)

      tst = b'*(c.*c.*c.*c) - (1/5);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 5th order condition A (tst = %g)\n', tst);
         end
      end

      tst = (b.*c.*c)'*(A*c) - (1/10);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 5th order condition B (tst = %g)\n', tst);
         end
      end

      tst = b'*((A*c).*(A*c)) - (1/20);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 5th order condition C (tst = %g)\n', tst);
         end
      end
      
      tst = (b.*c)'*A*(c.*c) - (1/15);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 5th order condition D (tst = %g)\n', tst);
         end
      end

      tst = b'*A*(c.*c.*c) - (1/20);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 5th order condition E (tst = %g)\n', tst);
         end
      end

      tst = (b.*c)'*A*A*c - (1/30);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 5th order condition F (tst = %g)\n', tst);
         end
      end

      tst = -(1/40);
      for i=1:s
         for j=1:s
            for k=1:s
               tst = tst + b(i)*A(i,j)*c(j)*A(j,k)*c(k);
            end
         end
      end
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 5th order condition G (tst = %g)\n', tst);
         end
      end
      
      tst = b'*A*A*(c.*c) - (1/60);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 5th order condition H (tst = %g)\n', tst);
         end
      end
   
   end

   tst = b'*A*A*A*c - (1/120);
   if (abs(tst) > tol)
      Ofail = true;
      Lfail = true;
      if (reportL>1)
         fprintf('    Method fails linear 5th order condition I (tst = %g)\n', tst);
      end
   end

   if (reportL>1)
      if (~Ofail),  fprintf('  Method passes order 5 conditions\n'); end
      if (~Lfail && Ofail),  fprintf('  Method passes linear order 5 condition\n'); end
   end
   if (~Ofail)
      q = 5;
   end
   if (~Lfail)
      lq = 5;
   end
end

% check for sixth order
if (~Ofail || ~Lfail) 

   if (~Ofail)

      tst = b'*(c.*c.*c.*c.*c) - (1/6);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition A (tst = %g)\n', tst);
         end
      end
      
      tst = (b.*c.*c.*c)'*A*c - (1/12);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition B (tst = %g)\n', tst);
         end
      end

      tst = -(1/24);
      for i=1:s
         for j=1:s
            for k=1:s
               tst = tst + b(i)*c(i)*A(i,j)*c(j)*A(j,k)*c(k);
            end
         end
      end
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition C (tst = %g)\n', tst);
         end
      end

      tst = (b.*c.*c)'*A*(c.*c) - (1/18);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition D (tst = %g)\n', tst);
         end
      end

      tst = -(1/36);
      for i=1:s
         for j=1:s
            for k=1:s
               tst = tst + b(i)*A(i,j)*c(j)*c(j)*A(i,k)*c(k);
            end
         end
      end
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition E (tst = %g)\n', tst);
         end
      end

      tst = (b.*c)'*A*(c.*c.*c) - (1/24);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition F (tst = %g)\n', tst);
         end
      end

      tst = b'*A*(c.*c.*c.*c) - (1/30);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition G (tst = %g)\n', tst);
         end
      end

      tst = (b.*c.*c)'*A*A*c - (1/36);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition H (tst = %g)\n', tst);
         end
      end

      tst = -(1/72);
      for i=1:s
         for j=1:s
            for k=1:s
               for l=1:s
                  tst = tst + b(i)*A(i,j)*A(j,k)*c(k)*A(i,l)*c(l);
               end
            end
         end
      end
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition I (tst = %g)\n', tst);
         end
      end

      tst = -(1/48);
      for i=1:s
         for j=1:s
            for k=1:s
               tst = tst + b(i)*c(i)*A(i,j)*c(j)*A(j,k)*c(k);
            end
         end
      end
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition J (tst = %g)\n', tst);
         end
      end

      tst = -(1/60);
      for i=1:s
         for j=1:s
            for k=1:s
               tst = tst + b(i)*A(i,j)*c(j)*c(j)*A(j,k)*c(k);
            end
         end
      end
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition K (tst = %g)\n', tst);
         end
      end

      tst = -(1/120);
      for i=1:s
         for j=1:s
            for k=1:s
               for l=1:s
                  tst = tst + b(i)*A(i,j)*A(j,k)*c(k)*A(j,l)*c(l);
               end
            end
         end
      end
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition L (tst = %g)\n', tst);
         end
      end

      tst = (b.*c)'*A*A*(c.*c) - (1/72);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition M (tst = %g)\n', tst);
         end
      end

      tst = -(1/90);
      for i=1:s
         for j=1:s
            for k=1:s
               tst = tst + b(i)*A(i,j)*c(j)*A(j,k)*c(k)*c(k);
            end
         end
      end
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition N (tst = %g)\n', tst);
         end
      end

      tst = b'*A*A*(c.*c.*c) - (1/120);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition O (tst = %g)\n', tst);
         end
      end

      tst = (b.*c)'*A*A*A*c - (1/144);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition P (tst = %g)\n', tst);
         end
      end

      tst = -(1/180);
      for i=1:s
         for j=1:s
            for k=1:s
               for l=1:s
                  tst = tst + b(i)*A(i,j)*c(j)*A(j,k)*A(k,l)*c(l);
               end
            end
         end
      end
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition Q (tst = %g)\n', tst);
         end
      end

      tst = -(1/240);
      for i=1:s
         for j=1:s
            for k=1:s
               for l=1:s
                  tst = tst + b(i)*A(i,j)*A(j,k)*c(k)*A(k,l)*c(l);
               end
            end
         end
      end
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition R (tst = %g)\n', tst);
         end
      end

      tst = b'*A*A*A*(c.*c) - (1/360);
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 6th order condition S (tst = %g)\n', tst);
         end
      end

   end

   tst = b'*A*A*A*A*c - (1/720);
   if (abs(tst) > tol)
      Ofail = true;
      Lfail = true;
      if (reportL>1)
         fprintf('    Method fails linear 6th order condition (tst = %g)\n', tst);
      end
   end

   if (reportL>1)
      if (~Ofail),  fprintf('  Method passes order 6 conditions\n'); end
      if (~Lfail && Ofail),  fprintf('  Method passes linear order 6 condition\n'); end
   end
   if (~Ofail)
      q = 6;
   end
   if (~Lfail)
      lq = 6;
   end
   
end


% check for linear seventh order
if (~Lfail)
   tst = b'*A*A*A*A*A*c - (1/factorial(7));
   if (abs(tst) > tol)
      Lfail = true;
      if (reportL>1)
         fprintf('    Method fails linear 7th order condition (tst = %g)\n', tst);
      end
   end

   if (reportL>1)
      if (~Lfail),  fprintf('  Method passes linear order 7 condition\n'); end
   end
   if (~Lfail)
      lq = 7;
   end
end

% check for linear eighth order
if (~Lfail)
   tst = b'*A*A*A*A*A*A*c - (1/factorial(8));
   if (abs(tst) > tol)
      Lfail = true;
      if (reportL>1)
         fprintf('    Method fails linear 8th order condition (tst = %g)\n', tst);
      end
   end

   if (reportL>1)
      if (~Lfail),  fprintf('  Method passes linear order 8 condition\n'); end
   end
   if (~Lfail)
      lq = 8;
   end
end

% report
if (reportL>0)
   fprintf('  Overall results:\n');
   fprintf('    order = %i,  linear order = %i\n', q,lq);
end



% check Kinnmark & Gray stability conditions
if (s >= 5)
   if (abs(sum(b)-1) < tol),              s1='T'; else, s1='F'; end
   if (abs(b'*c-1/2) < tol),              s2='T'; else, s2='F'; end
   if (abs(b'*A*c-1/6) < tol),            s3='T'; else, s3='F'; end
   if (abs(b'*A*A*c-1/(2*15)) < tol),     s4='T'; else, s4='F'; end
   if (abs(b'*A*A*A*c-1/(5*2*15)) < tol), s5='T'; else, s5='F'; end
   if (reportL>0)
      fprintf('    5-stage Kinnmark & Gray stability condition:');
      fprintf('      %s %s %s %s %s\n', s1, s2, s3, s4, s5)
   end
end
if (s >= 6)
   if (abs(sum(b)-1) < tol),                 s1='T'; else, s1='F'; end
   if (abs(b'*c-1/2) < tol),                 s2='T'; else, s2='F'; end
   if (abs(b'*A*c-1/6) < tol),               s3='T'; else, s3='F'; end
   if (abs(b'*A*A*c-1/24) < tol),            s4='T'; else, s4='F'; end
   if (abs(b'*A*A*A*c-2/(15*24)) < tol),     s5='T'; else, s5='F'; end
   if (abs(b'*A*A*A*A*c-2/(6*15*24)) < tol), s6='T'; else, s6='F'; end
   if (reportL>0)
      fprintf('    6-stage Kinnmark & Gray stability condition:');
      fprintf('      %s %s %s %s %s %s\n', s1, s2, s3, s4, s5, s6)
   end
end
if (s >= 7)
   if (abs(sum(b)-1) < tol),                     s1='T'; else, s1='F'; end
   if (abs(b'*c-1/2) < tol),                     s2='T'; else, s2='F'; end
   if (abs(b'*A*c-1/6) < tol),                   s3='T'; else, s3='F'; end
   if (abs(b'*A*A*c-4/(3*35)) < tol),            s4='T'; else, s4='F'; end
   if (abs(b'*A*A*A*c-4/(5*3*35)) < tol),        s5='T'; else, s5='F'; end
   if (abs(b'*A*A*A*A*c-8/(9*35*35)) < tol),     s6='T'; else, s6='F'; end
   if (abs(b'*A*A*A*A*A*c-8/(7*9*35*35)) < tol), s7='T'; else, s7='F'; end
   if (reportL>0)
      fprintf('    7-stage Kinnmark & Gray stability condition:'); 
      fprintf('      %s %s %s %s %s %s %s\n', s1, s2, s3, s4, s5, s6, s7)
   end
end
if (s >= 8)
   if (abs(sum(b)-1) < tol),                        s1='T'; else, s1='F'; end
   if (abs(b'*c-1/2) < tol),                        s2='T'; else, s2='F'; end
   if (abs(b'*A*c-1/6) < tol),                      s3='T'; else, s3='F'; end
   if (abs(b'*A*A*c-1/24) < tol),                   s4='T'; else, s4='F'; end
   if (abs(b'*A*A*A*c-1/(3*48)) < tol),             s5='T'; else, s5='F'; end
   if (abs(b'*A*A*A*A*c-1/(6*3*48)) < tol),         s6='T'; else, s6='F'; end
   if (abs(b'*A*A*A*A*A*c-4/(21*48*48)) < tol),     s7='T'; else, s7='F'; end
   if (abs(b'*A*A*A*A*A*A*c-4/(8*21*48*48)) < tol), s8='T'; else, s8='F'; end
   if (reportL>0)
      fprintf('    8-stage Kinnmark & Gray stability condition:');
      fprintf('      %s %s %s %s %s %s %s %s\n', s1, s2, s3, s4, s5, s6, s7, s8)
   end
end


% generate plot of stability region
if (doPlot)
   figure()
   xl = box(1:2);  yl = box(3:4);
   xax = plot(linspace(xl(1),xl(2),10),zeros(1,10),'k:'); hold on
   yax = plot(zeros(1,10),linspace(yl(1),yl(2),10),'k:');
   [X,Y] = stab_region2(A,b,box);
   plot(X,Y,'r-')
   set(get(get(xax,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
   set(get(get(yax,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
   axis(box)
   xlabel('Re(z)')
   ylabel('Im(z)')
   title(sprintf('%s stability region, order %i',mname,q))
   print(sprintf('%s_stab_region.png', fname), '-dpng');
   print(sprintf('%s_stab_region.eps', fname), '-depsc');
   savefig(sprintf('%s_stab_region.fig', fname));
end

% end of function
