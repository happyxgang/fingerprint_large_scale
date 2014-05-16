function [R,L] = match_query(D,SR,IX)
% [R,L] = match_query(D,SR,IX)
%     Match landmarks from an audio query against the database.
%     Rows of R are potential maxes, in format
%      songID  modalDTcount modalDT
%     i.e. there were <modalDTcount> occurrences of hashes 
%     that occurred in the query and reference with a difference of 
%     <modalDT> frames.
%     Input D is either audio waveform at sample rate SR, or a
%     3-column array of hashes (e.g., from landmark2hash, or get_hashes_for_track).
%     L returns the actual landmarks that this implies for IX'th return.
% 2008-12-29 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3
  IX = 1;
end


if size(D, 2) == 3
  % We were passed hashes, not audio
  Lq = [];
  Hq = D;

else

  % Target query landmark density
  % (reference is 7 lm/s, query can be maybe 4x denser?)
  dens = 20;

  % collapse stereo
  if size(D,2) == 2
    D = mean(D,2);
  end
  global xzg
    if xzg > 0
        disp(['using xzg find landmark']);
        %Rt = get_hash_hits(landmark2hash(find_landmarks(D,SR)));
        Lq = xzg_find_landmarks(D,SR, dens);
        %%Lq = fuzzify_landmarks(Lq);
        % Augment with landmarks calculated half-a-window advanced too
        landmarks_hopt = 0.032;
        %Lq = [Lq;xzg_find_landmarks(D(round(landmarks_hopt/4*SR):end),SR, dens)];
        %Lq = [Lq;xzg_find_landmarks(D(round(landmarks_hopt/2*SR):end),SR, dens)];
        Lq = [Lq;xzg_find_landmarks(D(round(3*landmarks_hopt/4*SR):end),SR, dens)];
        % add in quarter-hop offsets too for even better recall
    else
        %Rt = get_hash_hits(landmark2hash(find_landmarks(D,SR)));
        Lq = find_landmarks(D,SR, dens);
        %%Lq = fuzzify_landmarks(Lq);
        % Augment with landmarks calculated half-a-window advanced too
        landmarks_hopt = 0.032;
        Lq = [Lq;find_landmarks(D(round(landmarks_hopt/4*SR):end),SR, dens)];
        Lq = [Lq;find_landmarks(D(round(landmarks_hopt/2*SR):end),SR, dens)];
        Lq = [Lq;find_landmarks(D(round(3*landmarks_hopt/4*SR):end),SR, dens)];
    end

  % add in quarter-hop offsets too for even better recall
  Hq = unique(landmark2hash(Lq), 'rows');
end

disp(['landmarks ',num2str(size(Lq,1)),' -> ', num2str(size(Hq,1)),' hashes']);
% Rt row ={songid,delta_t,hash}
Rt = get_hash_hits(Hq);
nr = size(Rt,1);

if nr > 0
  % Find all the unique tracks referenced
  [utrks,xx] = unique(sort(Rt(:,1)),'first');
  utrkcounts = diff([xx',nr]);

  [utcvv,utcxx] = sort(utrkcounts, 'descend');
  % Keep at most 20 per hit
  %  utcxx = utcxx(1:min(20,length(utcxx)));
  
  % 最大的track的hash  count
  utrkcounts = utrkcounts(utcxx);
  % 最大的track id
  utrks = utrks(utcxx);
  
  nutrks = length(utrks);
  R = zeros(nutrks,4);
  % 对于每一个最多的hash的track
  for i = 1:nutrks
    % 找出所有的songid delta time 序列
    tkR = Rt(Rt(:,1)==utrks(i),:);
    % Find the most popular time offset
    [dts,xx] = unique(sort(tkR(:,2)),'first');
    dtcounts = 1+diff([xx',size(tkR,1)]);
    [vv,xx] = max(dtcounts);
    %R(i) = [songid, 最大delta_t_num, delta_t,song_id匹配的hash数目]
    R(i,:) = [utrks(i),vv,dts(xx),size(tkR,1)];
    R(i,:) = [utrks(i),sum(abs(tkR(:,2)-dts(xx))<=1),dts(xx),size(tkR,1)];
  end

  % Sort by descending match count
  [vv,xx] = sort(R(:,2),'descend');
  R = R(xx,:);

  % Extract the actual landmarks
  H = Rt((Rt(:,1)==R(IX,1)) & (Rt(:,2)==R(IX,3)),:);
  % Restore the original times
  for i = 1:size(H,1)
    hix = find(Hq(:,3)==H(i,3));
    hix = hix(1);  % if more than one...
    H(i,2) = H(i,2)+Hq(hix,2);
    L(i,:) = hash2landmark(H(i,:));
  end


  % Return no more than 100 hits, and only down to 10% the #hits in
  % most popular
  maxrtns = 100;
  if size(R,1) > maxrtns
    R = R(1:maxrtns,:);
  end
  maxhits = R(1,2);
  nuffhits = R(:,2)>(0.1*maxhits);
  %R = R(nuffhits,:);

else
  R = [];
  disp('*** NO HITS FOUND ***');
end
