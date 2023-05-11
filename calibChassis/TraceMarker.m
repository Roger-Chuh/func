function [curMarkerId, curMarkerPt] = TraceMarker(prvMarkerId, prvMarkerPt, curMarkerIdCandidate, curMarkerPtCandidate, prvImgL, curImgL, threshold)
%% input

% prvMarkerId          : n x 1 vector, stands for n different marker ids in previous frame.

% prvMarkerPt          : n x 8 matrix, stands for n sets of marker points detected in previous frame,
%                        each row stands for: [x0 y0 x1 y1 x2 y2 x3 y3]

% curMarkerIdCandidate : n x 1 vector, stands for n different marker ids in
%                        current frame, for those whose ids aren't identified, assign '-1' to
%                        corresponding entry.
%                        for example: 
%                                  curMarkerIdCandi = [1; 2; -1; 4; 5]; means marker 1 2 4 5
%                                  are successfully identified, but there are also ONE marker's corner
%                                  points are detected without an id assigned to it.
% 

% curMarkerPtCandidate : n x 8 matrix, stands for n sets of marker points detected in current frame,
%                        each row stands for: [x0 y0 x1 y1 x2 y2 x3 y3]

% prvImgL              : previous left image

% curImgL              : current left image

% threshold            : a threshold(unit: pixel) to determine whether it is an valid id.





%% output

% curMarkerId          : n x 1 vector, stands for n different marker ids in current frame.

% curMarkerPt          : n x 8 matrix, stands for n sets of marker points detected in current frame,
%                        each row stands for: [x0 y0 x1 y1 x2 y2 x3 y3]


prvPtX = prvMarkerPt(:,1:2:end);
prvPtY = prvMarkerPt(:,2:2:end);

curPtCandiX = curMarkerPtCandidate(:,1:2:end);
curPtCandiY = curMarkerPtCandidate(:,2:2:end);

[lkPtList, inTrackFlag] = LK(prvImgL, curImgL, [prvPtX(:) prvPtY(:)]);
lkPtList(~inTrackFlag, :) = -1;

curPtLkX = reshape(lkPtList(:,1), [], 4);
curPtLkY = reshape(lkPtList(:,2), [], 4);

alignedIdPair = [];
validId = [];
for i = 1 : size(curMarkerPtCandidate, 1)
    
    curId = curMarkerIdCandidate(i);
    
    if curId == -1
        % current corners without id
        error = 1000*ones(size(prvMarkerPt, 1), 1);
        for j = 1 : size(prvMarkerPt, 1)
            goodLK = find(curPtLkX(j,:) > 0);
            error(j,1) = mean(sqrt((curPtCandiX(i, goodLK) - curPtLkX(j, goodLK)).^2 + (curPtCandiY(i, goodLK) - curPtLkY(j, goodLK)).^2));
        end
        [minError,  minId] = min(error);
        if minError < threshold
            alignedIdPair = [alignedIdPair; [i prvMarkerId(minId)]];
        end
    else
        % current corners with id
        id = find(prvMarkerId == curId);
        if ~isempty(id)
            % old id reDetected
             goodLK = find(curPtLkX(id,:) > 0);
            err = mean(sqrt((curPtCandiX(i, goodLK) - curPtLkX(id, goodLK)).^2 + (curPtCandiY(i, goodLK) - curPtLkY(id, goodLK)).^2));
            if err < threshold
                validId = [validId; i];
            end
        else
            % new id detected
            validId = [validId; i];
        end
    end
end


curMarkerId = [alignedIdPair(:,2); curMarkerIdCandidate(validId, :)];
curMapingrkerPt = [curMarkerIdCandidate(alignedIdPair(:,1), :); curMarkerIdCandidate(validId, :)];

end