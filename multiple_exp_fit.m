function [y] = multiple_exp_fit(x,tmax,lcc,s,varargin)
    % Connects multiple experiments
    % x,tmax,lcc,s;c,l1,ld,lg,k2n,k2p,k3p,ksn,Ks,kg,a,A,st,rt,numOfExp

    c = varargin{1};
    l1 = varargin{2};
    ld = varargin{3};
    lg = varargin{4};
    k2n = varargin{5};
    k2p = varargin{6};
    k3p = varargin{7};
    ksn = varargin{8};
    Ks = varargin{9};
    kg = varargin{10};
    a = varargin{11};
    A = varargin{12};
    st = varargin{13};
    rt = varargin{14};
    numOfExp = varargin{15};

    beginidx = 1;
    xlen = size(x,1);

    y = zeros(xlen,1);

    for i=1:numOfExp
        endidx = selectTime(x,tmax*i);
        
        if  x(beginidx) >= tmax*(i-1) && x(beginidx) < tmax*i
            if x(endidx) < tmax*i 
                y(beginidx:endidx) = exp_fit(x(beginidx:endidx)-tmax*(i-1),st,rt,l1,ld,lg,k2n,k2p,k3p,ksn,Ks,0,kg,s(:,i),a,lcc(:,i),[],7,c,A,0);
            else
                y(beginidx:endidx-1) = exp_fit(x(beginidx:endidx-1)-tmax*(i-1),st,rt,l1,ld,lg,k2n,k2p,k3p,ksn,Ks,0,kg,s(:,i),a,lcc(:,i),[],7,c,A,0);
            end
        end
        beginidx = endidx;
    end

end
