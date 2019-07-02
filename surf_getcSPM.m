function cSPM=surf_getcSPM(type,varargin)
% function MAP=surf_getcSPM(type,varargin)
% generates a statistical map from columns of functional gifti files 
% INPUT: 
%   type: Type of linear model 
%       'onesample_t':  tests the Hypothesis that the mean effect is bigger
%                       than zero
%                       the only implemented type of map right now is the
%       'regression':   simple or multiple regression  
%                       'regression',X 
%                       'regression',X,'no_intercept'
%                       where X is a regression matrix 
%                       The intercept term will be added, unless otherwise
%                       specified
%   varargin: 
%       'data',matrix: 
%                       Specify data in a matrix format. 
%       'Z_P':          adds to each contrast a Normal value, corresponding to the
%                       p-value of the contrast
%       'delta':        add to each contrast the effect-size (mean/SD) for
%                       comparision of effects with different number of
%                       subjects 
%        'maskthreshold',p: only computes the statistics when n>=p*N subject are 
%                       measured at this node 
%                       set maskthreshold to 1 if you want to compute the
%                       t-value only if all subjects are measured 
%        'contrasts',C: Contrast vectors (one contrast per row) 
%                       if not specified, assume identity matrix 
% OUTPUT: 
%       cSPM.

% _________________________________________________________________
% Adapted from caret_getcSPM function (EBerlot, May 2019)
% jdiedric@uwo.ca
if(nargin<1 || isempty(type)) 
    type='onesample_t';
end; 
C=[]; 
cSPM.data=[];
cSPM.N=0;
switch (type)
    case 'onesample_t'
        c=1;
    case 'regression'
        X=varargin{1};
        c=2;
    otherwise 
        error(['Unknown type: ' type]);
end;
maskthreshold=0.7; % Compute only when at least 70% of subjects are measured 
z_p=0;
delta=0;
is_intercept=1;
while c<=length(varargin)
    switch (varargin{c})
        case 'no_intercept'
            is_intercept=0;    
            c=c+1;
        case 'data'   
            cSPM.data=varargin{c+1};
            c=c+2;
        case 'Z_P'             % Compute and write out the z-value corresponding to the P-value
            z_p=1;
            c=c+1;
        case 'delta'            % Compute and write out effect-size
            delta=1;
            c=c+1;
        case 'maskthreshold'
            maskthreshold=varargin{c+1};
            c=c+2;
        case 'contrasts' 
            C=varargin{c+1}; 
            c=c+2; 
        otherwise 
            error(['Unknown option:' varargin{c}]);
    end;
end;
N=size(cSPM.data,2);  % Number of available data points 

% Check how much data I have at each Node to be able to use the
% maskthreshold 
% This assumes that 0 and NaNs are missing values 
miss=sum(isnan(cSPM.data) | cSPM.data==0.0,2);
cSPM.N=N-miss; 
indx=find(isnan(cSPM.data));
cSPM.data(indx)=0;
switch (type)
    case 'onesample_t'
        X=ones(N,1);cSPM.X=X;
        cSPM.b=nanmean(cSPM.data')';
        cSPM.ResVar=(nanstd(cSPM.data')').^2;
        cSPM.ResVar(cSPM.ResVar==0)=NaN;
        cSPM.bvar=inv(X'*X);
        cSPM.P=1;
        cSPM.title='one-sample t-test: mean>0';
    case 'regression'
        if (is_intercept)
            X=[ones(N,1) X];
        end;
        cSPM.X=X;
        iX=inv(X'*X);
        cSPM.b=(iX*X'*cSPM.data')';
        res=(cSPM.data'-cSPM.X*cSPM.b')';
        cSPM.bvar=iX;
        cSPM.P=rank(X);
        cSPM.ResVar=sum(res.^2,2)./(cSPM.N-cSPM.P);
        clear res; 
        cSPM.title='regression t-test: b>0';
end;

% Make contrasts 
if (isempty(C)) 
    C=eye(size(X,2));
end;
for i=1:size(C,1)
    cSPM.con(i).df=[1 max(cSPM.N)-cSPM.P];
    cSPM.con(i).STAT='T';
    c=C(i,:)'; 
    cSPM.con(i).con=c'*cSPM.b'; 
    cSPM.con(i).Z=(c'*cSPM.b'./sqrt((c'*cSPM.bvar*c)*cSPM.ResVar'))';  
    if (maskthreshold>0)
        cSPM.con(i).Z(cSPM.N<maskthreshold*N)=NaN;
        cSPM.con(i).con(cSPM.N<maskthreshold*N)=NaN;
    end;
    cSPM.con(i).Z(isnan(cSPM.con(i).Z))=0;
    if (z_p)
        cSPM.con(i).Z_P=norminv(tcdf(cSPM.con(i).Z,cSPM.con(i).df(2)));
    end;
    if (delta)
        cSPM.con(i).delta=c'*cSPM.b(:,i)./sqrt(cSPM.ResVar);
    end;
end;

