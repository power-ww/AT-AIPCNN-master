function varargout = ndwt2_working(X,level,varargin)
%NDWT2 Nondecimated 2-D wavelet transform.
%   NDWT2 will be removed in a future release of MATLAB. Use the
%   following function instead:
%       <a href="matlab:help swt2">swt2</a>

% Error in R2015a
% error(message('Wavelet:warnobsolete:ErrorReplaceNDWT2'));
nbIn = length(varargin);
if nbIn < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
elseif nbIn > 5
    error(message('MATLAB:narginchk:tooManyInputs'));
end

LoD = cell(1,2); HiD = cell(1,2); LoR = cell(1,2); HiR = cell(1,2);
if ischar(varargin{1})
    [LD,HD,LR,HR] = wfilters(varargin{1}); 
    for k = 1:2
        LoD{k} = LD; HiD{k} = HD; LoR{k} = LR; HiR{k} = HR;
    end

elseif isstruct(varargin{1})
    if isfield(varargin{1},'w1') && isfield(varargin{1},'w2')
        for k = 1:2
            [LoD{k},HiD{k},LoR{k},HiR{k}] = ...
                wfilters(varargin{1}.(['w' int2str(k)]));
        end
    elseif isfield(varargin{1},'LoD') && isfield(varargin{1},'HiD') && ...
           isfield(varargin{1},'LoR') && isfield(varargin{1},'HiR')
        for k = 1:2
            LoD{k} = varargin{1}.LoD{k}; HiD{k} = varargin{1}.HiD{k};
            LoR{k} = varargin{1}.LoR{k}; HiR{k} = varargin{1}.HiR{k};
        end
    else
        error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'));
    end
        
elseif iscell(varargin{1})
    if ischar(varargin{1}{1})
        for k = 1:2
            [LoD{k},HiD{k},LoR{k},HiR{k}] = wfilters(varargin{1}{k});
        end
    else
        LoD(1:end) = varargin{1}(1); HiD(1:end) = varargin{1}(2);
        LoR(1:end) = varargin{1}(3); HiR(1:end) = varargin{1}(4);
    end
else
    
end
nextArg = 2;

dwtEXTM = 'sym';
while nbIn>=nextArg
    argName = varargin{nextArg};
    argVal  = varargin{nextArg+1};
    nextArg = nextArg + 2;
    switch argName
        case 'mode' , dwtEXTM = argVal;
    end
end

% Initialization.
if isempty(X) , varargout{1} = []; return; end
sX = size(X);
X = double(X);
sizes = zeros(level+1,length(sX));
sizes(level+1,:) = sX;

for k=1:level
    dec = decFUNC(X,LoD,HiD,dwtEXTM);
    X = dec{1,1,1};
    sizes(level+1-k,:) = size(X);
    dec = reshape(dec,4,1,1);
    if k>1
        cfs(1) = [];
        cfs = cat(1,dec,cfs);
    else
        cfs = dec;
    end
end

WT.sizeINI = sX;
WT.level = level;
WT.filters.LoD = LoD;
WT.filters.HiD = HiD;
WT.filters.LoR = LoR;
WT.filters.HiR = HiR;
WT.mode = dwtEXTM;
WT.dec = cfs;
WT.sizes = sizes;
varargout{1} = WT;