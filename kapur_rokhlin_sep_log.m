function c = kapur_rokhlin_sep_log(k)
% Compute Kapur-Rokhlin correction weights for the "separable" logarithmic 
% singularity of the form 
%       -\log|x| * \phi(x), x\in[-a,a]
% Use Table 7 of [1] when weights available for the given k,
% otherwise compute weights via the zeta correction approach [2].
%
% Note that the output weight is negative of the K-R weights, and the
% center weight c(1) is scaled by 2.
%
% [1] Kapur, S., & Rokhlin, V. (1997). High-order corrected trapezoidal 
%     quadrature rules for singular functions. SIAM J. on Numer. Anal., 
%     34(4), 1331-1356.
% [2] Wu, B., & Martinsson, P.G. (2020). Zeta Correction: A New Approach 
%     to Constructing Corrected Trapezoidal Rule for Singular Integral 
%     Operators. (arxiv)
%     
% Bowei Wu, May 2020; Updated Jan 2021.

if ismember(k,[1,2,3,4,5,8,16,20,21]) % extract precomputed weights
    c = precom_weights(k); 
else
    % zeta connection approach (require Symbolic Math Toolbox!)
    rhs = imag(zeta(-2*vpa(0:(k-1))'+1i*eps)/eps); % derivative of zeta ia complex step differentiation
    V = vpa(0:(k-1)).^((0:2:2*(k-1))'); % extended precision Vandermonde matrix to combat ill-conditioning for large k
    c = double(V\rhs);
end
% post-processing
c = -c'; % negative of K-R weights; row vector
c(1) = 2*c(1); % scale center weight by 2


function c = precom_weights(k)
% precomputed weights
switch k
    case 1
        c = -.9189385332046728E+00;
    case 2
        c = [-.8884900761462795E+00;
             -.3044845705839327E-01];
    case 3
        c = [-.8788820090191140E+00;
             -.4325921322794724E-01;
             0.3202689042388492E-02];
    case 4
        c = [-.8742261022238798E+00;
             -.5024307342079858E-01;
             0.5996233119529027E-02;
             -.4655906795234226E-03];
    case 5
        c = [-.8714860605482545E+00;
             -.5462714010179898E-01;
             0.8188266460029230E-02;
             -.1091885919666338E-02
             0.7828690501786438E-04];
    case 8
        c = [-.8674568230266233E+00;
             -.6148248184914681E-01;
             0.1238115647253341E-01;
             -.2897732763288696E-02;
             0.6052303442088134E-03;
             -.9784299004952013E-04;
             0.1051436637368180E-04;
             -.5537586803550720E-06];
    case 16
        c = [-.8641774186256508E+00;
             -.6745551530557688E-01;
             0.1689299430945688E-01;
             -.5726331418330602E-02;
             0.2079973982393914E-02;
             -.7402722529894703E-03;
             0.2463303649153491E-03;
             -.7432197238379932E-04;
             0.1984260034353227E-04;
             -.4580836910292848E-05;
             0.8916453665775555E-06;
             -.1418448246845963E-06;
             0.1767037741742959E-07;
             -.1614096203394717E-08;
             0.9602551109604562E-10;
             -.2789306492547372E-11];
    case 20
        c = [-.8635304116647681E+00;
             -.6867878881201404E-01;
             0.1792611880692437E-01;
             -.6505172477564357E-02;
             0.2603300521146428E-02;
             -.1053029105125390E-02;
             0.4121099047922528E-03;
             -.1519803191376791E-03;
             0.5184904980367620E-04;
             -.1612303155465844E-04;
             0.4509136652480968E-05;
             -.1119040467513331E-05;
             0.2428371655709533E-06;
             -.4528845255185268E-07;
             0.7103184896000958E-08;
             -.9102854426840433E-09;
             0.9146251630822981E-10;
             -.6753249440965090E-11;
             0.3256755515337396E-12;
             -.7693371141612778E-14];
    case 21
        c = [-.8634075135173825E+00
             -.6891288052131983E-01;
             0.1812828891950665E-01;
             -.6663392565672225E-02;
             0.2715373083556168E-02;
             -.1124755545067623E-02
             0.4534905432204645E-03
             -.1734369464708259E-03
             0.6181105535120864E-04
             -.2024524074674083E-04
             0.6020613356244510E-05
             -.1606613597759635E-05
             0.3799671084527263E-06
             -.7853207506864614E-07
             0.1394746012004667E-07
             -.2083589766806165E-08
             0.2544214502140801E-09
             -.2437043148484079E-10
             0.1716505712892348E-11
             -.7901799480102858E-13
             0.1783115591485395E-14];
end