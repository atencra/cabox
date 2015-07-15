function fx = ringach_malone_func(a,x)
% ringach_malone_func   Nonlinearity function from Ringach and Malone (2007)
%
% fx = ringach_malone_func(a,x)
% --------------------------------------------------------
% x : independent variable, a vector of projection values
%
% a(1) : A, or gain
% a(2) : theta, or threshold
% a(3) : sigman, or external noise
%
% fx : function evaluated at x for parameters a.
%
% The nonlinearity function is:
%
% fx(z|A,theta,sigman) = [] * [] + *** exp()
%
% The function is defined in Ringach and Malone (2007), page 7675. Miller and
% Troyer (2002) and Hansel and van Vreeswijk (2002) derived the function 
% theoretically.
%
% To fit this function to a nonlinearity, make the following call:
%
% fit(x, y, fittype('ringach_malone_func(a,x)'));
%
% Here, x is the abscissa values, y is the estimated nonlinearity values, 
% and fittype specifies the function object that is to be fitted.
%
% caa 3/8/10
%
% fx = ringach_malone_func(a,x)

%{
A = a(1)
theta = a(2)
sigma = a(3)

part1 = [' (A .* ( x-theta ) ./ 2) .* (1 + myerf( (x-theta)./sigma./sqrt(2), 0)) +'];
part3 = ['(A .* sigma ./ sqrt(2*pi) ) .* exp( -( x - theta ).^2 ./ sqrt(2) ./ sigma.^2 )'];
ringach_func = [part1 part2 part3 part4];
%}


fx = (a(1) .* ( x-a(2) ) ./ 2) .* (1 + erfcore( (x-a(2))./a(3)./sqrt(2), 0)) + ...
(a(1) .* a(3) ./ sqrt(2*pi) ) .* exp( -( x - a(2) ).^2 ./ sqrt(2) ./ a(3).^2 );

% The Ringach and Malone (2007) paper doesn't have the baseline parameter, a(4):
%fx = (a(1) .* ( x-a(2) ) ./ 2) .* (1 + myerf( (x-a(2))./a(3)./sqrt(2), 0)) + ...
%(a(1) .* a(3) ./ sqrt(2*pi) ) .* exp( -( x - a(2) ).^2 ./ sqrt(2) ./ a(3).^2 );

return;







