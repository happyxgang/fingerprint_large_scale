function H = get_hashes_for_track(ID)
% H = get_hashes_for_track(ID)
%     Return all the hashes held in the global hash table for the
%     specified track ID
% 2013-04-24 Dan Ellis dpwe@ee.columbia.edu

global HashTable HashTableCounts
[nhtcols,nhtrows] = size(HashTable);

TIMESIZE=16384;

matches = find(((HashTable(:)-TIMESIZE/2)/TIMESIZE) == ID);

nh = length(matches);
H = zeros(nh, 3);

for i = 1:length(matches)
  hash = floor((matches(i)-1)/nhtcols);
  time = rem(HashTable(matches(i)), TIMESIZE);
  H(i,:) = [ID, time, hash];
end

[vv,ix] = sort(H(:,2));
H = H(ix,:);
