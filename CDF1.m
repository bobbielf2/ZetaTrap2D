function cd = CDF1(k)
% central difference coefficients for the 1st derivative
% Input:
%   k: integer, such that 2k+1 points about the center are used, order of
%      accuracy is 2k
% Output:
%   cd: central difference coefficients at nodes 1:k
%       (coefficients at nodes -k:-1 are antisymmetric, coeff at 0 is 0)

if k <= 22
    cd = precom_coeffs(k);
else
    x = sym(1:k);               % half of the -n:n grid
    A=2*x.^((1:2:2*k)');        % odd moments Vandermonde, times 2 bc symmetry
    b = zeros(k,1); b(1) = 1;   % rhs: 1st deriv of monomials eval at 0
    cd = double(A\b);          	% sym central diff coeffs at +/-1, +/-2, ...
end
cd = cd'; % row vector

function c = precom_coeffs(k)
% precomputed coeffs
switch k
    case 0
        c = 0;
    case 1
        c = 0.5;
    case 2
        c = [8;-1]/12;
    case 3
        c = [45;-9;1]/60;
    case 4
        c = [0.8
            -0.2
             4/105
            -1/280];
    case 5
        c = [ 8.333333333333334e-01
            -2.380952380952381e-01
             5.952380952380952e-02
            -9.920634920634920e-03
             7.936507936507937e-04];
    case 6
        c = [8.571428571428571e-01
            -2.678571428571428e-01
             7.936507936507936e-02
            -1.785714285714286e-02
             2.597402597402597e-03
            -1.803751803751804e-04];
    case 7
        c = [8.750000000000000e-01
            -2.916666666666667e-01
             9.722222222222222e-02
            -2.651515151515152e-02
             5.303030303030303e-03
            -6.798756798756799e-04
             4.162504162504163e-05];
    case 8
        c = [8.888888888888888e-01
            -3.111111111111111e-01
             1.131313131313131e-01
            -3.535353535353535e-02
             8.702408702408702e-03
            -1.554001554001554e-03
             1.776001776001776e-04
            -9.712509712509713e-06];
    case 9
        c = [9.000000000000000e-01
            -3.272727272727273e-01
             1.272727272727273e-01
            -4.405594405594405e-02
             1.258741258741259e-02
            -2.797202797202797e-03
             4.495504495504496e-04
            -4.627725215960510e-05
             2.285296402943462e-06];
    case 10
        c = [9.090909090909091e-01
            -3.409090909090909e-01
             1.398601398601399e-01
            -5.244755244755245e-02
             1.678321678321678e-02
            -4.370629370629371e-03
             8.814714697067638e-04
            -1.285479226655697e-04
             1.202787580496559e-05
            -5.412544112234514e-07];
    case 11
        c = [9.166666666666666e-01
            -3.525641025641026e-01
             1.510989010989011e-01
            -6.043956043956044e-02
             2.115384615384616e-02
            -6.221719457013574e-03
             1.481361775479422e-03
            -2.728824323251568e-04
             3.638432431002091e-05
            -3.118656369430363e-06
             1.288700979103456e-07];
    case 12
        c = [9.230769230769231e-01
            -3.626373626373626e-01
             1.611721611721612e-01
            -6.799450549450549e-02
             2.559793148028442e-02
            -8.295625942684766e-03
             2.245432585989862e-03
            -4.911883781852822e-04
             8.316416985147635e-05
            -1.020651175449937e-05
             8.068388738734680e-07
            -3.081676254377829e-08];
    case 13
        c = [9.285714285714286e-01
            -3.714285714285714e-01
             1.702380952380952e-01
            -7.510504201680672e-02
             3.004201680672269e-02
            -1.054105852867463e-02
             3.162317558602388e-03
            -7.905793896505971e-04
             1.597130080102216e-04
            -2.499855777551295e-05
             2.840745201762835e-06
            -2.083213147959412e-07
             7.396023010506791e-09];
    case 14
        c = [9.333333333333333e-01
            -3.791666666666667e-01
             1.784313725490196e-01
            -8.178104575163399e-02
             3.443412452700378e-02
            -1.291279669762642e-02
             4.216423411469851e-03
            -1.173890608875129e-03
             2.722065180000299e-04
            -5.103872212500561e-05
             7.423814127273543e-06
            -7.852111096154709e-07
             5.368964852071596e-08
            -1.780524058084968e-09];
    case 15
        c = [9.375000000000000e-01
            -3.860294117647059e-01
             1.858660130718954e-01
            -8.804179566563468e-02
             3.873839009287926e-02
            -1.537237702098383e-02
             5.390314020344980e-03
            -1.640530354018037e-03
             4.253226843750468e-04
            -9.186969982501010e-05
             1.606113633304372e-05
            -2.181141971154086e-06
             2.157173378064481e-07
            -1.381441079548682e-08
             4.297816691929233e-10];
    case 16
        c = [9.411764705882353e-01
            -3.921568627450980e-01
             1.926384588923289e-01
            -9.391124871001032e-02
             4.293085655314757e-02
            -1.788785689714482e-02
             6.666282073470121e-03
            -2.187373805357383e-03
             6.221863268572113e-04
            -1.507605330461704e-04
             3.045667334266069e-05
            -4.985467362637911e-06
             6.347544652695483e-07
            -5.894148606074377e-08
             3.549164752044786e-09
            -1.039794360950621e-10];
    case 17
        c = [9.444444444444444e-01
            -3.976608187134503e-01
             1.988304093567251e-01
            -9.941520467836257e-02
             4.699627857522595e-02
            -2.043316459792432e-02
             8.027314663470271e-03
            -2.809560132214594e-03
             8.644800406814137e-04
            -2.305280108483770e-04
             5.239272973826750e-05
            -9.936552191740388e-06
             1.528700337190829e-06
            -1.831622523823113e-07
             1.602669708345224e-08
            -9.106077888325135e-10
             2.520713602304536e-11];
    case 18
        c = [9.473684210526315e-01
            -4.026315789473684e-01
             2.045112781954887e-01
            -1.045796308954204e-01
             5.092573330559601e-02
            -2.298731017266487e-02
             9.457636185324973e-03
            -3.501144164759725e-03
             1.152640054241885e-03
            -3.334423014056882e-04
             8.362189376945157e-05
            -1.788579394513270e-05
             3.195476833869862e-06
            -4.636294513427255e-07
             5.245100863675278e-08
            -4.338778287966682e-09
             2.333460591847627e-10
            -6.121733034168158e-12];
    case 19
        c = [9.500000000000000e-01
            -4.071428571428571e-01
             2.097402597402597e-01
            -1.094297007340486e-01
             5.471485036702428e-02
            -2.553359683794466e-02
             1.094297007340486e-02
            -4.255599472990777e-03
             1.486082355647573e-03
            -4.611979724423503e-04
             1.257812652115501e-04
            -2.975470789950647e-05
             6.008162172015729e-06
            -1.014365042028630e-06
             1.392265743960864e-07
            -1.491713297100926e-08
             1.169971213412491e-09
            -5.972826014418121e-11
             1.489070197500363e-12];
    case 20
        c = [9.523809523809523e-01
            -4.112554112554113e-01
             2.145680406549972e-01
            -1.139892715979672e-01
             5.836250705815924e-02
            -2.805889762411502e-02
             1.247062116627334e-02
            -5.066189848798544e-03
             1.863426151282223e-03
            -6.149306299231337e-04
             1.803315630273119e-04
            -4.649173109297885e-05
             1.040374402080646e-05
            -1.988951062801234e-06
             3.182321700481975e-07
            -4.143648047502572e-08
             4.216112480765733e-09
            -3.143592639167433e-10
             1.527251484615757e-11
            -3.627222275962422e-13];
    case 21
        c = [9.545454545454546e-01
            -4.150197628458498e-01
             2.190382081686429e-01
            -1.182806324110672e-01
             6.186986926117361e-02
            -3.055302185736968e-02
             1.402944881205751e-02
            -5.926232687851878e-03
             2.282697035320724e-03
            -7.952621929504456e-04
             2.485194352970143e-04
            -6.903317647139285e-05
             1.686783497491047e-05
            -3.580111913042222e-06
             6.497240138484033e-07
            -9.877561021343969e-08
             1.223227371064269e-08
            -1.184892610147725e-09
             8.418973808944358e-11
            -3.901475667559581e-12
             8.846883599908347e-14];
    case 22
        c = [9.565217391304348e-01
            -4.184782608695652e-01
             2.231884057971014e-01
            -1.223244147157191e-01
             6.523968784838351e-02
            -3.300817539947975e-02
             1.560977752881801e-02
            -6.829277668857879e-03
             2.741502146638288e-03
            -1.002361722364624e-03
             3.313592470626857e-04
            -9.827075709457099e-05
             2.591756231065609e-05
            -6.016576964973735e-06
             1.214156072210916e-06
            -2.096815585232667e-07
             3.036113064590287e-08
            -3.584300145696867e-09
             3.312831970348837e-10
            -2.247993122736711e-11
             9.957887586873580e-13
            -2.160285530210178e-14];
end