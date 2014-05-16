function [N,T] = add_tracks(D,SR,ID)
% [N,T] = add_tracks(D,SR,ID)
%    Add one or more tracks to the hashtable database.
%    <D, SR> define the waveform of the track, and ID is its
%    reference ID.
%    If D is a char array, load that wavefile.  Second arg is ID.
%    If D is a cell array, load each of those wavefiles; second arg
%    is vector of IDs.
%    N returns the total number of hashes added, T returns total
%    duration in secs of tracks added.
% 2008-12-29 Dan Ellis dpwe@ee.columbia.edu

% Target query landmark density
% (reference is 7 lm/s)
dens = 20;

if isnumeric(D)
  H = landmark2hash(find_landmarks(D,SR,dens),ID);
%  save_hashes(H);
  record_hashes(H);
  N = length(H);
  T = length(D)/SR;
elseif ischar(D)
  if nargin < 3
    ID = SR;
  end
  [D,SR] = wavread(D);
  [N,T] = add_tracks(D,SR,ID);
elseif iscell(D)
  
  disp(['Target density = ',num2str(dens),' hashes/sec']);
  
  nd = length(D);
  if nargin < 3
    if nargin < 2
      % omitting IDs defaults to track number
      ID = 1:nd;
    else
      ID = SR;
    end
  end
  N = 0;
  T = 0;
  for i = 1:nd
    disp(['Adding #',num2str(ID(i)),' ',D{i},' ...']);
    [n,t] = add_tracks(D{i},ID(i));
    N = N + n;
    T = T + t;
  end
  disp(['added ',num2str(nd),' tracks (',num2str(T),' secs, ', ...
        num2str(N),' hashes, ',num2str(N/T),' hashes/sec)']);
else
  error('I cant tell what D is');
end

%
% 2010-04-20  Building index for artist20 on hog @ 20 hash/sec:
% Adding #1413 /proj/hog3/labrosa/data/artistid/artist20/mp3s-32k/u2/Zooropa/10-The_Wanderer.mp3 ...
% added 1413 tracks (348315.3467 secs, 7182012 hashes, 20.6193 hashes/sec)
% Elapsed time is 8667.995757 seconds.
% >> save HTA20-20hps HashTable HashTableCounts Names
% -rw-r--r-- 1 dpwe 1512 31812971 Apr 20 22:18 HTA20-20hps.mat
%
% with dens = 10
%added 1413 tracks (348315.3467 secs, 3051605 hashes, 8.761 hashes/sec)
%Elapsed time is 8230.792979 seconds.
