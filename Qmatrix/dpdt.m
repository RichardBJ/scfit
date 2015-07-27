function dpdt = dpdt(t, p, tf, wfm, Q, C)
%DPDT Returns the derivative of the state probabilities with time for an
%ion channel mechanism described by Q matrix
%   Use DPDT with an ode solver to give the channel state probabilities in
%   response to a concentration waveform.
%
%   t - time at which to give the derivatives of p
%   p - vector of state probabilities at time t
%   tf - time vector corresponding to the waveform
%   wfm - concentration waveform
%   Q - Q matrix describing the mechanism
%   C - matrix same size as Q identifying concentration-dependent states
%       0 indicates no concentration-dependence and a positive integer is
%       the index of the ligand whose concentration the rate depends on

    wfm = interp1q (tf, wfm, t);
%     dpdt = zeros(size(Q,1),1);

    conc = ones(size(C));
    conc(C~=0) = wfm(C(C~=0));
    conc = prod(conc,3);
    
    rates = conc.*Q;
    rates = rates';
    rates = rates+diag(-sum(rates));
    dpdt = rates*p;
end