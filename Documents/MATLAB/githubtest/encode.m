function [ lpary, plocs, flocs ] = lowpass_encode( grayimary, log2N, ovlp, cfq )
%LOWPASS_ENCODE Image low-pass transform encoding.
    
    % PART 1: Zeropadding and getting patch locations - we do this for you
    %===========================================================
    N = 2^(log2N);                  % get patch sizes
    [im_expd, plocs] = expdimg( grayimary, N, ovlp );
    % plocs is the location of the left top of each patch in the original
    % grayimage. 
    numps = size(plocs,1); % number of patches
    
	% PART 2: Take the DFT and throw away high freq. components
	%===========================================================
    % identify frequency locations to be retained
	
	% Q - Determine the the zero frequency location after fftshift:
    %%%% FILL IN THE VALUE BELOW %%%%
    
    % zero frequency location is the location (round(n/2) , round(n/2)), as
    % described above (for each patch)
    % since plocs is left top location for each patch, we need to add it by
    % N/2
	zeroidx = plocs(:) + N/2;
    zeroidx = reshape(zeroidx, size(plocs));
	
	% Use meshgrid to get the row and column indices of the N x N grid
    %   row indices first, then column indices: MATLAB order.
	[fidxs(:,:,2), fidxs(:,:,1)] = meshgrid(1:N, 1:N);
    
	% Q - Determine the frequency locations to be kept:
	%   i)  Produce a logical array fkeep of the locations satisfying the
	%       requirements from the handout. Use fkeep to determine numfs.
    %%%% FILL IN THE CODE BELOW %%%%

	% we need to calculate the region for each patch
    % actually we could only do it for the first patch, then apply the
    % result to the other patches
    first_center = zeroidx(1,:);
    fkeep = zeros(N,N);
    for i=1:N
        for j=1:N
            if (i-first_center(1))^2 + (j-first_center(2))^2 < cfq * floor(N/2)^2
                fkeep(i,j) = 1;
            end
        end
    end
    numfs = sum(sum(fkeep));		% number of frequency locations kept
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %   ii) Use fidxs and fkeep to produce flocs.
    %%%% FILL IN THE CODE BELOW %%%%
    
    flocs = zeros(numfs, 2);
    cnt = 1;
    for i = 1:size(fkeep,2)
        for j = 1:size(fkeep,1)
            if fkeep(j,i) ~= 0
                flocs(cnt,:) = [j, i];
                cnt = cnt + 1;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
	% PART 3: For each patch, produce a col. of retained freq. values
	%===========================================================
    lpary = zeros(numfs,numps);
    flidxs = sub2ind([N, N], flocs(:,1), flocs(:,2));
    for i = 1:numps
        % Q - Extract the patch given the patch location:
        %%%% FILL IN THE VALUE BELOW %%%%
        patch = im_expd(plocs(i,1):plocs(i,1)+N-1, plocs(i,2):plocs(i,2)+N-1);
        
        % Q - For the fftshift(fft2()) of each patch, produce a column 
        % of retained frequency values using flocs or fkeep:
        %%%% FILL IN THE CODE BELOW %%%%
        
        fft_patch = fftshift(fft2(patch));
        temp = fft_patch(flidxs);
        lpary(:,i) = temp;
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end



%% DO NOT MODIFY THIS SUBFUNCTION:
%===============================================================
function [ im_expd, plocs ] = expdimg(grayimary, N, ovlp)
    M = N - 2*ovlp;

    % enlargement for subpatch partition:
    imsize = size(grayimary);
    Mpartsize = ceil(imsize/M)*(M);      
    
    % zeropad, including overlap:
    im_expd = zeros(Mpartsize  + 2*ovlp);
    im_expd(ovlp+(1:imsize(1)), ovlp+(1:imsize(2))) = grayimary;
    
    % get the locations of the patches:
    [plocs(:,:,2), plocs(:,:,1)] = meshgrid(1:M:Mpartsize(2), 1:M:Mpartsize(1));
    plocs = reshape(plocs, numel(plocs(:,:,1)), 2);
end