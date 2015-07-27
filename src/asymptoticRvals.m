function [ s, areaR, r, c, Wprime, mu, a ] = asymptoticRvals( q, A, F, tau ) %#codegen
%ASYMPTOTICRVALS Calculate the values used for the asymptotic approximation
%of the exact pdf of dwell times with missed events 
%   See Hawkes et al (1992) Phil Trans R Soc Lond B pp.393-394
%
% The asymptotic behavior of the matrix function, R, depends on values of
% s that render the matrix W(s) singular. The values of s that render
% W(s) singular are the roots of the determinant of W(s). That is,
%       det W(s) = 0
% Jalali and Hawkes (1992) Adv Appl Prob 24 pp.302-321 prove that if the
% Q matrix is irreducible and reversible (which should hold if the gating
% mechanism follows microscopic reversibility - see Colquhoun and Hawkes
% 1982 pp. 24-25), then det W(s) = 0 has exactly kA roots

tol=1e-9;
kA = length(A);
eqFFt = expm(q(F,F)*tau);
uf = ones(length(F),1);
s = inf(kA,1);
r = zeros(kA,kA);
c = zeros(kA,kA);
mu = -1./s;
Wprime = zeros(kA,kA,kA);
areaR = zeros(kA,kA,kA);
a = zeros(kA,1);
npts = 100*(kA) + 1;
x0 = zeros(1,npts);
y0 = zeros(1,npts);
idx1 = zeros(1,kA);
idx2 = zeros(1,kA);

%Find the kA roots of det W(s) = 0
%First, guess that the roots are the same as the tau^-1 for the uncorrected
%distribution
% the real() function ensures proper compiling with Matlab Coder
tempx = real(sort([0;eig(q(A,A))]));

% Start with the values of tempx closest to zero, which correspond to time
% constants with the longest duration since tau ~= 1 ./ s
% This will allow us to (hopefully) find the time constants that are above
% the imposed resolution
% KKO 27 Aug 2014
tempx = flipud(tempx);

% What if all eigenvalues are the same??? Does det W still have kA roots?
% Is the Q matrix then reducible?

% Construct initial guesses by evenly spacing 100 points in between the
% values in tempx
% Then find when the sign of detW(x0) changes and save that index for later
% input to rootsdetW
uu = 1;
for ii=1:npts
    j1 = ceil(ii/100);
    j2 = mod(ii-1,100);
    if ii==npts
        x0(ii) = tempx(j1);
    else
        x0(ii) = (tempx(j1+1)-tempx(j1)).*(j2./100) + tempx(j1);
    end
    y0(ii) = real(detW(x0(ii),q,A,F,tau));
    if ii>1
        if y0(ii)>0 && y0(ii-1)<0
            idx1(uu) = ii-1;
            idx2(uu) = ii;
            uu = uu+1;
        elseif y0(ii)<0 && y0(ii-1)>0
            idx1(uu) = ii-1;
            idx2(uu) = ii;
            uu = uu+1;
        end
    end
    if uu>kA
        break
    end
end

ii=1;
x0 = real(x0);
for rr=1:uu-1
    idx = [idx1(rr), idx2(rr)];
    tmp = rootsdetW(q,A,F,tau,x0(idx));
    if all(abs(tmp-s)>tol)
        s(ii)=tmp;
        ii=ii+1;
    end
end

% If not all roots have been found, try using fzeroKO (in rootsdetW) using
% the first and last values of x0 as inputs
if any(isinf(s))
    ii = find(isinf(s),1);
    tmp = rootsdetW(q,A,F,tau,x0(end));
    if all(abs(tmp-s)>tol)
        s(ii)=tmp;
    end
end
if any(isinf(s))
    ii = find(isinf(s),1);
    tmp = rootsdetW(q,A,F,tau,x0(1));
    if all(abs(tmp-s)>tol)
        s(ii)=tmp;
    end
end

% If we still haven't found all roots, then give up
if any(isinf(s))
    return
end

% Calculate the right and left eigenvectors of H(s), where s are the roots
% of det W(s) = 0, which are also eigenvalues of H(s)
%
% The left eigenvector, r, is a solution to rW(s) = 0, ru=1, where u is a
% vector of ones - this is similar to finding equilibrium vector of a Q
% matrix (e.g. Hawkes and Sykes 1990)
%
% The right eigenvector, c, is a solution to c'W(s)' = 0, c'u =1

for ii=1:kA
    r(ii,:) = [zeros(1,kA), 1]/[W(s(ii),q,A,F,tau), ones(kA,1)];
    c(ii,:) = [zeros(1,kA), 1]/[W(s(ii),q,A,F,tau)', ones(kA,1)];
end
c = c';

mu = -1./s;
for ii=1:kA
    Wprime(:,:,ii) = dWds(s(ii),q,A,F,tau);
    areaR(:,:,ii) = c(:,ii)*r(ii,:) ./ (r(ii,:)*Wprime(:,:,ii)*c(:,ii));
    a(ii) = mu(ii)*phi(q,A,F,tau)*c(:,ii)*r(ii,:)*q(A,F)*eqFFt*uf ./ (r(ii,:)*Wprime(:,:,ii)*c(:,ii));
end

end

