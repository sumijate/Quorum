function ret = exp_fit(varargin)
    % scale the data
    %   t, st, rt, l1,  ld,  lg,  k2n, k2p, k3p, ksn,  Ks, bn,  kg,    s,    a, divtime,              IC, ret,    c,   A,  tau
    %   1,  2,  3,  4,   5,   6,    7,   8,   9,  10,  11, 12,  13,   14,   15,      16,              17,  18,   19,  20,   21
    % 0:1, 10,  1, 20, 0.5, 0.7, 1000, 0.5, 100, 0.1, 400,  0, 1.5, 1000, 2300,      [], [0,0,0,0,0,0,0], 1:7,    0, 500, 0.25

    Defaults = {0:1,10,1,20,0.5,0.7,1000,0.5,100,0.1,400,0,1.5,1000,2300,[],[0,0,0,0,0,0,0],1:7,0,500,0.25};
    idx = ~cellfun('isempty',varargin);
    Defaults(idx) = varargin(idx);
    
    t = Defaults{1};
    tau = Defaults{21};
        
    dt = 64/(768*5);
    tall = (0):dt:(max(t)+tau);
            
    % Parameters
    st = Defaults{2};
    rt = Defaults{3};
    l1 = Defaults{4}; % 1/h
    ld = Defaults{5};
    lg = Defaults{6};
    k2n = Defaults{7};
    k2p = Defaults{8};
    k3p = Defaults{9};
    ksn = Defaults{10};
    Ks = Defaults{11};
    bn = Defaults{12};
    kg = Defaults{13}; % 1/h
    s = Defaults{14};
    a = Defaults{15};
    DivTime = Defaults{16};
    IC = Defaults{17};
    selRet = Defaults{18};
    c = Defaults{19}; % 1/h
    A = Defaults{20};

    X = exp_ode_solver(tall,st,rt,l1,ld,lg,k2n,k2p,k3p,ksn,Ks,bn,kg,s,a,DivTime,IC,1:7);

    ret = A*X(:,selRet)+c;

    t = t + tau;
    ret = ret(selectTime((tall)',t'),:);
end