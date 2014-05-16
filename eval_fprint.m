function [S,R,QF] = eval_fprint(Q,SR,T,DUR,SNR)
% [S,R,QF] = eval_fprint(Q,SR,T,DUR,SNR)
%    Evaluate the fingerprinting system over a set of queries.
%    Q is a cell array of query waveforms, each at sampling rate
%    SR.  T is the ground-truth track indices that should be 
%    returned (0 => not found).  Return S as the proportion of
%    queries correctly identified.  R is a matrix of actual 
%    top-hit results, with 4 columns: track_id nmatch t_offs total_match,
%    as returned by match_query.
%    QF returns the actual truncated & noised queries from Q.
%    DUR truncates all queries to this many seconds (default all).
%    SNR adds noise at this SNR (in dB) relative to query (default 60).
% 2010-04-21 DAn Ellis dpwe@ee.columbia.edu

if nargin < 4;  DUR = 999; end
if nargin < 5;  SNR = 60; end

nq = length(Q);
s = 0;

maxdursamps = round(DUR*SR);

for i = 1:nq
  dd = Q{i};
  if length(dd) > maxdursamps
    dd = dd(1:maxdursamps);
  end
  % figure noise level given actual energy
  noise = ( (10^(-SNR/20))*sqrt(mean(dd.^2)) )*randn(length(dd),1);

  dd = dd + noise;

  r = match_query(dd,SR);
  R(i,:) = r(1,:);
  QF{i} = dd;
end

if nargin > 2
  S = mean(R(:,1)==T');
else
  S = 0;
end


  