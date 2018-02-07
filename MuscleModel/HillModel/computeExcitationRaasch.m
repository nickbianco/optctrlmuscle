function e = computeExcitationRaasch(a, vA, tauDeact, tauAct)
% Computes muscle excitation from muscle activation a and derivative of
% activation vA.

% See De Groote et al. (2009) equation (17) (Raasch)

td = (ones(size(a,1),1)*tauDeact);
ta = (ones(size(a,1),1)*tauAct);

e = zeros(size(a));
idx = vA<=0;
e(idx) = td(idx) .* vA(idx) + a(idx);

c1 = 1./ta - 1./td;
c2 = 1./td;
D = (c2 + c1 .* a).^2 + 4*c1.*vA;
idx = vA>0;
e(idx) = (a(idx) .* c1(idx) - c2(idx) + sqrt(D(idx)))./(2*c1(idx));

end

