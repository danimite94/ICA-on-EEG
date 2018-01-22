function [g, mu_opt] = kurt_gradient_optstep(f, X, s, P)

% Computes optimal step size in the gradient-based optimization of the normalized kurtosis contrast
% (single iteration).
%
% Data-based version.
%
% See reference below for details.
%
%
% SYNTAX: [g, mu_opt] = kurt_gradient_optstep(f, X, s, P); 
%
%
% OUTPUTS:
%         g      : search direction (gradient vector)
%
%         mu_opt : optimal step size globally optimizing the normalized kurtosis contrast function
%                  along direction g from f.
%
%
% INPUTS:
%         f      : current extracting vector coefficients
%
%         X      : sensor-output data matrix (one signal per row, one sample per column)
%
%         s      : source kurtosis sign; if zero, the maximum absolute value of the contrast is sought
%
%         P      : projection matrix (used in deflationary orthogonalization; identity matrix otherwise).
%
%
% REFERENCE:
%
% - V. Zarzoso, P. Comon and M. Kallel,  <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/eusipco-2006.pdf.gz">"How Fast is FastICA?"</a>, 
%   in: Proceedings EUSIPCO-2006, XIV European Signal Processing Conference, 
%   Florence, Italy, September 4-8, 2006. 
%
%
% Please, report any bugs, comments or suggestions to <a href = "mailto:zarzoso@i3s.unice.fr">zarzoso(a)i3s.unice.fr</a>.
%
%
% HISTORY:
% 
%    <Please add modification date here>: - <please add modification details here>
%
%    2008/03/26: - added this help
%
%    2008/03/25: - projecting the gradient after normalization seems to improve conditioning
%                  and accelerate convergence in the extraction of the last sources
%
%    2008/03/24: - created by Vicente Zarzoso.



[L, T] = size(X);

mu_opt = 0;  % default optimal step-size value


%%% Compute search direction (gradient vector)

% compute necessary interim values

y = f'*X;

ya2 = y.*conj(y);
y2 = y.*y;
ya4 = ya2.*ya2;

Eya2 = mean(ya2);
Ey2 = mean(y2);
Eya4 = mean(ya4);


if abs(Eya2) < eps  % check for zero denominator
    g = zeros(L, 1);
    
else % compute gradient if contrast denominator is not null
    Eycx = X*y'/T;
    Eyx = X*y.'/T;
    Ey3x = X*(ya2.*y)'/T;

    % contrast numerator and denominator at current point
    p1 = Eya4 - 2*Eya2^2 - abs(Ey2)^2;
    p2 = Eya2;

    g = 4*( (Ey3x - 2*Eya2*Eycx - Eyx*Ey2') - Eycx*p1/p2 )/p2^2;
       
    if norm(g) > eps
        g = g/norm(g);      % normalize the gradient -> parameter of interest: direction 
    end                     % improves conditioning of opt step-size polynomial  

    g = P*g;                % project if required (after normalization)
     
    
%%% Compute optimal step size

gg = g'*X;

% calculate interim values

ya2 = y.*conj(y);
ga2 = gg.*conj(gg);
reygc = real(y.*conj(gg));

g2 = gg.*gg;
yg = y.*gg;

Eya2reygc = mean(ya2.*reygc);
Ereygc2 = mean(reygc.^2);
Ega2reygc = mean(ga2.*reygc);
Ega4 = mean(ga2.^2);
Eya2ga2 = mean(ya2.*ga2);

Ega2 = mean(ga2);
Ereygc = mean(reygc);

Eg2 = mean(g2);
Eyg = mean(yg);


% obtain contrast-function polynomials

p11 = [Ega4, 4*Ega2reygc, 4*Ereygc2+2*Eya2ga2, 4*Eya2reygc, Eya4];

p2 = [Ega2, 2*Ereygc, Eya2];    % square-root of denominator
p12 = conv(p2, p2);             % denominator

p13 = [Eg2, 2*Eyg, Ey2];
p13 = conv(p13, conj(p13));

p1 = p11 - 2*p12 - p13;         % numerator

Jnum = p1;
Jden = p12;


% compute derivatives

p1d = [4, 3, 2, 1].*p1(1:4);
p2d = [2, 1].*p2(1:2);

% derivative's numerator
Jd_num = conv(p1d, p2) - 2*conv(p2d, p1);   % 5th-degree coefficient, Jd_num(1), is always 0
Jd_den = conv(p2, p12);
 
if (Jd_num*Jd_num')/(Jd_den*Jd_den') > eps  % check if contrast is constant (e.g., because only one source is left;
                                            % in that case, any extracting vector is valid up to scale)
rr = roots(Jd_num);

n = polyval(Jnum, rr);
d = polyval(Jden, rr);

nonzero_d = find(abs(d) > eps);     % check roots not shared by denominator

if ~isempty(nonzero_d)     
    
   n = n(nonzero_d);
   d = d(nonzero_d);
   rr = rr(nonzero_d);
 
   Jkm_val = n./d;                  % normalized kurtosis
   
   if s
    Jkm_val = real(s*Jkm_val);      % maximize or minimize kurtosis value, depending on kurtosis sign
   else
    Jkm_val = abs(Jkm_val);         % maximize absolute kurtosis value, if no sign is given
   end

   [Jmax, im] = max(Jkm_val);   
                                      
   mu_opt = real(rr(im));           % keep only real part of the first root, in the case more than one are found

end % if ~isempty(nonzero_d)
end % if (Jd_num*Jd_num')/(Jd_den*Jd_den')

end % if abs(Eya2) < eps



return;




%%% Alternative equivalent expressions of optimal step size polynomial (see EUSIPCO-2006 paper)

a = abs(y).^2; b = abs(g).^2; c = real(y.*conj(g));    
d = y.^2;      e = g.^2;      f = y.*g;

am = mean(a); a2m = mean(a.^2);
bm = mean(b); b2m = mean(b.^2); abm = mean(a.*b);
cm = mean(c); c2m = mean(c.^2); acm = mean(a.*c); bcm = mean(b.*c);
dm = mean(d); em = mean(e);     fm = mean(f);

h0 = a2m - 2*am^2 - abs(dm)^2;
h1 = 4*acm - 8*am*cm - 4*real(dm*fm');
h2 = 4*c2m + 2*abm - 8*cm^2 - 4*am*bm - 4*abs(fm)^2 - 2*real(dm*em');
h3 = 4*bcm - 8*bm*cm - 4*real(em*fm');
h4 = b2m - 2*bm^2 - abs(em)^2;

p1_alt = [h4, h3, h2, h1, h0]

i0 = am; i1 = 2*cm; i2 = bm;

p2_alt = [i2, i1, i0]

mu0 = -2*h0*i1 + h1*i0;
mu1 = -4*h0*i2 - h1*i1 + 2*h2*i0;
mu2 = -3*h1*i2 + 3*h3*i0;
mu3 = -2*h2*i2 + h3*i1 + 4*h4*i0;
mu4 = -h3*i2 + 2*h4*i1;

Jdnum_alt = [0, mu4, mu3, mu2, mu1, mu0]

average_square_difference = sum((Jd_num - Jdnum_alt).^2)/6  








