function [root, fx, ea, iter] = falsePosition(func, xl, xu, es, maxit, varargin)
%falsePosition finds the root of a function using false position method
%Created by Tim Hunt
%Created on 9/30/19
%INPUTS:
%  func - the function to be analysed
%  xl - the lower bound for the guess of the root
%  xu - the upper bound for the guess of the root
%  es - the maximum estimated percent error of the root
%  maxit - the maximum number of iterations
%OUTPUTS:
%  root - the root calculated by the function
%  fx - the value of the function at the root
%  ea - the calculated maximum error of the root from the real root
%  iter - the number of iterations used

iter = 0;
lastR = xl;

if nargin < 5
    maxit = 200;
    if nargin < 4
        es = 0.0001;
    end
end

if func(xl)*func(xu) >= 0
    error('there is no sign change between func at xl and func at xu');
end

while 1
   iter = iter + 1;
   
   fl = func(xl);fu = func(xu);
   
   root = xu - (fu*(xu-xl))/(fu-fl);
   
   fx = func(root);
   
   if fx == 0
       ea = 0;
       break
   end
   
   ea = 100*abs((root-lastR)/root);
   
   if ea <= es || iter == maxit || fx == 0
       break
   elseif fx*fl > 0
       xl = root;
   else 
       xu = root;
   end
   lastR = root;
end

end