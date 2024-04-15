function X = exp_ode_solver(varargin)
    % Model for GFP
    %   t, st, rt, l1,  ld,  lg,  k2n, k2p, k3p, ksn,    Ks, bn,  kg,    s,    a, divtime,              IC, ret
    %   1,  2,  3,  4,   5,   6,    7,   8,   9,  10,    11, 12,  13,   14,   15,      16,              17,  18
    % 0:1, 10,  1, 20, 0.5, 0.7, 1000, 0.5, 100, 0.1,   400,  0, 1.5, 1000, 2300,      [], [0,0,0,0,0,0,0], 1:7
    
    Defaults = {0:1,10,1,20,0.5,0.7,1000,0.5,100,0.1,400,0,1.5,1000,2300,[],[0,0,0,0,0,0,0],1:7};
    idx = ~cellfun('isempty',varargin);
    Defaults(idx) = varargin(idx);

    %% Parameters
    t = Defaults{1};
    st = Defaults{2};
    rt = Defaults{3};
    l1 = Defaults{4}; % 1/h
    l2 = Defaults{5};
    l3 = l2;
    l4 = l2;
    lg = Defaults{6};
    ln = lg;
    k2n = Defaults{7};
    k2p = Defaults{8};
    k3n = k2n;
    k3p = Defaults{9};
    k4n = k3n;
    k4p = k3p;
    ksn = Defaults{10};
    Ks = ksn/Defaults{11};
    bn = Defaults{12};
    kg = Defaults{13}; % 1/h
    s = Defaults{14};
    a = Defaults{15};
    b = (1000 - a)/(log(2)/3);
    divtime = Defaults{16};
    IC = Defaults{17};
    ret = Defaults{18};
    
    %% Time span
    tinitial = t(1);
    tfinal = t(end);
    tspan = [tinitial, tfinal];

    %% Model/ system of differential equations
    Model = @(t, X) [(a+b*divtime(floor(t*12)+1))*rt + 2*k2n*X(2) - 2*k2p*X(1)*X(1) - (l1+divtime(floor(t*12)+1))*X(1);
                    k2p*X(1)*X(1) + k3n*X(3) - 2*k3p*s(floor(t*12)+1)*X(2) - (k2n+l2+divtime(floor(t*12)+1))*X(2);
                    2*k3p*s(floor(t*12)+1)*X(2) + 2*k4n*X(4) - k4p*s(floor(t*12)+1)*X(3) - (k3n+l3+divtime(floor(t*12)+1))*X(3);
                    k4p*s(floor(t*12)+1)*X(3) - (2*k4n+l4+divtime(floor(t*12)+1))*X(4);
                    Ks*X(4)*(st-X(5)) - (ksn+divtime(floor(t*12)+1))*X(5);
                    bn*st + ((a+b*divtime(floor(t*12)+1))-bn)*X(5) - (kg+ln+divtime(floor(t*12)+1))*X(6);
                    kg*X(6) - (lg+divtime(floor(t*12)+1))*X(7)];

    %% Solution
    sol = ode45(Model, tspan, IC);
    X = deval(sol, t);

    X = X(ret,:)';

end