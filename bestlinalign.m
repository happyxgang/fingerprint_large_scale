function [O,S] = bestlinalign(d1,sr1,d2,sr2)
% [O,S] = bestlinalign(d1,sr1,d2,sr2)
%    Use Shazam-style fingerprints to find best time alignment
%    between two waveforms d1 (at sample rate sr1) and d2 (at
%    sample rate sr2).  O returns the offset (in seconds) by which
%    d1 should be delayed to make it synchronize with d2; S returns
%    the time scaling (as a factor) which should additionally be
%    applied to d1 to make it line up over the entire duration; S
%    == 1.001 means that d2 is 0.1% slower than d1.
% 2010-04-27 Dan Ellis dpwe@ee.columbia.edu

timestep = 0.032;

if ischar(d1)
  [P1,N1] = fileparts(d1);
  [d1,sr1] = mp3read(d1,0,1);
else
  N1 = 'd1';
end

if ischar(d2)
  [P2,N2] = fileparts(d2);
  [d2,sr2] = mp3read(d2,0,1);
else
  N2 = 'd2';
end


dens = 10;
% make <id timestep hash> rows for both audio
H1 = double(landmark2hash(find_landmarks(d1,sr1,dens)));
H2 = double(landmark2hash(find_landmarks(d2,sr2,dens)));

nH1 = size(H1,1);
nH2 = size(H2,1);

% full match array: rows: d1 cols: d2
MA = (repmat(H1(:,3),1,nH2) == repmat(H2(:,3)',nH1,1));

% indices of matches - <index in H1> <index in H2>
MPs = [rem(find(MA)-1,nH1)+1,floor((find(MA)-1)/nH1)+1];

disp([num2str(length(MPs)),' common hashes found']);

% Calculate time differences (delay of d2 relative to d1)
DTs = H2(MPs(:,2),2) - H1(MPs(:,1),2);

% Time threshold
tthr = 0.2;
% Find modal time skew, with overlapping 400ms bins
minDT = min(DTs);
maxDT = max(DTs);
DTbins = minDT:tthr:maxDT;
[n,x] = hist(DTs,DTbins/timestep);
% Overlap adjacent bins
n(2:end) = n(2:end)+n(1:end-1);
% Find mode
[vv,xx] = max(n);
bestDT = DTbins(xx);

% % assume good matches are tightly distributed around median - allow
% % +/-200ms (or +/- 20 s from center for 1% time skew)
% vv = sort(DTs);
% bestDT = vv(round(length(vv)/2));
disp(['Best time skew = ',num2str(bestDT),' s']);
goodDTs = find( abs(DTs*timestep-bestDT) < tthr );
disp([num2str(length(goodDTs)),' matches found within ',num2str(tthr),' s of best']);

% Least-squares fit to those points
A = [H1(MPs(goodDTs,1),2),ones(length(goodDTs),1)] \ H2(MPs(goodDTs,2),2);
O = A(2) * timestep;
S = A(1);

disp(['Best match for time T sec in ',N1,' is ',...
      num2str(O),'+',num2str(S),'*T in ',N2]);

plot(H1(MPs(:,1),2)*timestep, timestep*DTs, '.b', ...
     H1(MPs(goodDTs,1),2)*timestep, timestep*DTs(goodDTs), '.r', ...
     H1(:,2)*timestep, O+(S-1)*H1(:,2)*timestep, '-r');
title(['Matching landmarks: ',N1,' vs. ',N2]);
xlabel(['time in ',N1,' / s']);
ylabel(['time in ',N2,' - time in ',N1,' / s']);