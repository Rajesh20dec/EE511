function Proj6Q3()
clc;
%clear all;
No_of_samples=1000;
Random_Sample=zeros(1,No_of_samples);

%--------below code is for alpha 0.05 and beta 0---------
Random_Sample=stblrnd(0.5,0,1,0,No_of_samples,1);
figure(1)
yyaxis left
xlim([-1500 1500]);
hist(Random_Sample,10000);
%h.Color = 'red';
title('Sample Vs Pdf(alpha=0.5)') 
xlabel('x-->')
ylabel('Frequency')
yyaxis right
t=-100:1:100;
p = stblpdf(t,0.5,0,1,0);
plot(t,p);
ylabel('PDF');
legend({'Hist of A','Theoretical PDF of alph stable dist.'},'FontSize',8)
figure(2)
plot([1:No_of_samples],Random_Sample)
title('Times series plot of smaples(Alpha=0.5,Beta=0.0)'); 
xlabel('sample Number')
ylabel('sample values')
sprintf('Range of  samples values collected at alha(0.5) is [%f %f]',min(Random_Sample),max(Random_Sample))

%--------below code is for alpha 1.0 and beta 0---------
figure(3)
Random_Sample=stblrnd(1,0,1,0,No_of_samples,1);
%hist(Random_Sample,20);
yyaxis left
xlim([-60 60]);
hist(Random_Sample,300);
%h.Color = 'red';
title('Sample Vs Pdf(alpha=1)') 
xlabel('x-->')
ylabel('Frequency')
yyaxis right
t=-50:0.1:50;
p = stblpdf(t,1,0,1,0);
plot(t,p);
ylabel('PDF');
legend({'Hist of A','Theoretical PDF of alph stable dist.'},'FontSize',8)

figure(4)
plot([1:No_of_samples],Random_Sample)
title('Times series plot of smaples(Alpha=1,Beta=0.0)');  
xlabel('sample Number')
ylabel('sample values')
sprintf('Range of  sample value collected at alha(1) is [%f %f]',min(Random_Sample),max(Random_Sample))

%--------below code is for alpha 1.8 and beta 0---------
figure(5)
Random_Sample=stblrnd(1.8,0,1,0,No_of_samples,1);
%hist(Random_Sample,20);
yyaxis left
xlim([-20 20]);
hist(Random_Sample,40);
%h.Color = 'red';
title('Sample Vs Pdf(alpha=1.8)') 
xlabel('x-->')
ylabel('Frequency')
yyaxis right
t=-50:0.1:50;
p = stblpdf(t,1.8,0,1,0);
plot(t,p);
ylabel('PDF');
legend({'Hist of A','Theoretical PDF of alph stable dist.'},'FontSize',8)

figure(6)
plot([1:No_of_samples],Random_Sample)
title('Times series plot of smaples(Alpha=1.8,Beta=0.0)');  
xlabel('sample Number')
ylabel('sample values')
sprintf('Range of  sample value collected at alha(1.8) is [%f %f]',min(Random_Sample),max(Random_Sample))


%--------below code is for alpha 2 and beta 0---------
figure(7)
Random_Sample=stblrnd(2,0,1,0,No_of_samples,1);
%hist(Random_Sample,20);
yyaxis left
xlim([-5 5]);
hist(Random_Sample,20);
%h.Color = 'red';
title('Sample Vs Pdf(alpha=2)') 
xlabel('x-->')
ylabel('Frequency')
yyaxis right
t=-50:0.1:50;
p = stblpdf(t,2,0,1,0);
plot(t,p);
ylabel('PDF');
legend({'Hist of A','Theoretical PDF of alph stable dist.'},'FontSize',8)

figure(8)
plot([1:No_of_samples],Random_Sample)
title('Times series plot of smaples(Alpha=2,Beta=0.0)');  
xlabel('sample Number')
ylabel('sample values')
sprintf('Range of  sample value collected at alha(2) is [%f %f]',min(Random_Sample),max(Random_Sample))



%--------below code is for alpha 0.5 and beta 0.75---------

figure(9)
yyaxis left

Random_Sample=stblrnd(0.5,0.75,1,0,No_of_samples,1);
xlim([-1500 1500]);
hist(Random_Sample,10000);
%h.Color = 'red';
title('Sample Vs Pdf(alpha=0.5 and Beta=0.75)') 
xlabel('x-->')
ylabel('Frequency')
yyaxis right
t=-100:0.1:100;
p = stblpdf(t,0.5,0.75,1,0);
plot(t,p);
ylabel('PDF');
legend({'Hist of A','Theoretical PDF of alph stable dist.'},'FontSize',8)

figure(10)
plot([1:No_of_samples],Random_Sample)
title('Times series plot of smaples(Alpha=0.5,Beta=0.75)'); 
xlabel('sample Number')
ylabel('sample values')

%--------below code is for alpha 1.0 and beta 0.75---------
figure(11)
Random_Sample=stblrnd(1,0.75,1,0,No_of_samples,1);
%hist(Random_Sample,20);
yyaxis left
xlim([-60 60]);
hist(Random_Sample,300);
%h.Color = 'red';
title('Sample Vs Pdf(alpha=1 and Beta=0.75)') 
xlabel('x-->')
ylabel('Frequency')
yyaxis right
t=-50:0.1:50;
p = stblpdf(t,1,0.75,1,0);
plot(t,p);
ylabel('PDF');
legend({'Hist of A','Theoretical PDF of alph stable dist.'},'FontSize',8)

figure(12)
plot([1:No_of_samples],Random_Sample)
title('Times series plot of smaples(Alpha=1,Beta=0.75)'); 
xlabel('sample Number')
ylabel('sample values')
%--------below code is for alpha 1.8 and beta 0.75---------
figure(13)
Random_Sample=stblrnd(1.8,0.75,1,0,No_of_samples,1);
%hist(Random_Sample,20);
yyaxis left
xlim([-20 20]);
hist(Random_Sample,40);
%h.Color = 'red';
title('Sample Vs Pdf(alpha=1.8 and Beta=0.75)') 
xlabel('x-->')
ylabel('Frequency')
yyaxis right
t=-50:0.1:50;
p = stblpdf(t,1.8,0.75,1,0);
plot(t,p);
ylabel('PDF');
legend({'Hist of A','Theoretical PDF of alph stable dist.'},'FontSize',8)

figure(14)
plot([1:No_of_samples],Random_Sample)
title('Times series plot of smaples(Alpha=1.8,Beta=0.75)'); 
xlabel('sample Number')
ylabel('sample values')
%--------below code is for alpha 2 and beta 0.75---------
figure(15)
Random_Sample=stblrnd(2,0.75,1,0,No_of_samples,1);
%hist(Random_Sample,20);
yyaxis left
xlim([-5 5]);
hist(Random_Sample,20);
%h.Color = 'red';
title('Sample Vs Pdf(alpha=2 and Beta=0.75)') 
xlabel('x-->')
ylabel('Frequency')
yyaxis right
t=-50:0.1:50;
p = stblpdf(t,2,0.75,1,0);
plot(t,p);
ylabel('PDF');
legend({'Hist of A','Theoretical PDF of alph stable dist.'},'FontSize',8)
figure(16)
plot([1:No_of_samples],Random_Sample)
title('Times series plot of smaples(Alpha=2,Beta=0.75)'); 
xlabel('sample Number')
ylabel('sample values')
end



function r = stblrnd(alpha,beta,gamma,delta,varargin)
%STBLRND alpha-stable random number generator.
% R = STBLRND(ALPHA,BETA,GAMMA,DELTA) draws a sample from the Levy 
% alpha-stable distribution with characteristic exponent ALPHA, 
% skewness BETA, scale parameter GAMMA and location parameter DELTA.
% ALPHA,BETA,GAMMA and DELTA must be scalars which fall in the following 
% ranges :
%    0 < ALPHA <= 2
%    -1 <= BETA <= 1  
%    0 < GAMMA < inf 
%    -inf < DELTA < inf
%
%
% R = STBLRND(ALPHA,BETA,GAMMA,DELTA,M,N,...) or 
% R = STBLRND(ALPHA,BETA,GAMMA,DELTA,[M,N,...]) returns an M-by-N-by-... 
% array.   
% 
%
% References:
% [1] J.M. Chambers, C.L. Mallows and B.W. Stuck (1976) 
%     "A Method for Simulating Stable Random Variables"  
%     JASA, Vol. 71, No. 354. pages 340-344  
%
% [2] Aleksander Weron and Rafal Weron (1995)
%     "Computer Simulation of Levy alpha-Stable Variables and Processes" 
%     Lec. Notes in Physics, 457, pages 379-392
%

if nargin < 4
    error('stats:stblrnd:TooFewInputs','Requires at least four input arguments.'); 
end

% Check parameters
if alpha <= 0 || alpha > 2 || ~isscalar(alpha)
    error('stats:stblrnd:BadInputs',' "alpha" must be a scalar which lies in the interval (0,2]');
end
if abs(beta) > 1 || ~isscalar(beta)
    error('stats:stblrnd:BadInputs',' "beta" must be a scalar which lies in the interval [-1,1]');
end
if gamma < 0 || ~isscalar(gamma)
    error('stats:stblrnd:BadInputs',' "gamma" must be a non-negative scalar');
end
if ~isscalar(delta)
    error('stats:stblrnd:BadInputs',' "delta" must be a scalar');
end


% Get output size
[err, sizeOut] = genOutsize(4,alpha,beta,gamma,delta,varargin{:});
if err > 0
    error('stats:stblrnd:InputSizeMismatch','Size information is inconsistent.');
end


%---Generate sample----

% See if parameters reduce to a special case, if so be quick, if not 
% perform general algorithm

if alpha == 2                  % Gaussian distribution 
    r = sqrt(2) * randn(sizeOut);

elseif alpha==1 && beta == 0   % Cauchy distribution
    r = tan( pi/2 * (2*rand(sizeOut) - 1) ); 

elseif alpha == .5 && abs(beta) == 1 % Levy distribution (a.k.a. Pearson V)
    r = beta ./ randn(sizeOut).^2;

elseif beta == 0                % Symmetric alpha-stable
    V = pi/2 * (2*rand(sizeOut) - 1); 
    W = -log(rand(sizeOut));          
    r = sin(alpha * V) ./ ( cos(V).^(1/alpha) ) .* ...
        ( cos( V.*(1-alpha) ) ./ W ).^( (1-alpha)/alpha ); 

elseif alpha ~= 1                % General case, alpha not 1
    V = pi/2 * (2*rand(sizeOut) - 1); 
    W = - log( rand(sizeOut) );       
    const = beta * tan(pi*alpha/2);
    B = atan( const );
    S = (1 + const * const).^(1/(2*alpha));
    r = S * sin( alpha*V + B ) ./ ( cos(V) ).^(1/alpha) .* ...
       ( cos( (1-alpha) * V - B ) ./ W ).^((1-alpha)/alpha);

else                             % General case, alpha = 1
    V = pi/2 * (2*rand(sizeOut) - 1); 
    W = - log( rand(sizeOut) );          
    piover2 = pi/2;
    sclshftV =  piover2 + beta * V ; 
    r = 1/piover2 * ( sclshftV .* tan(V) - ...
        beta * log( (piover2 * W .* cos(V) ) ./ sclshftV ) );      
          
end
    
% Scale and shift
if alpha ~= 1
   r = gamma * r + delta;
else
   r = gamma * r + (2/pi) * beta * gamma * log(gamma) + delta;  
end

end


%====  function to find output size ======%
function [err, commonSize, numElements] = genOutsize(nparams,varargin)
    try
        tmp = 0;
        for argnum = 1:nparams
            tmp = tmp + varargin{argnum};
        end
        if nargin > nparams+1
            tmp = tmp + zeros(varargin{nparams+1:end});
        end
        err = 0;
        commonSize = size(tmp);
        numElements = numel(tmp);

    catch
        err = 1;
        commonSize = [];
        numElements = 0;
    end
end


function p = stblpdf(x,alpha,beta,gam,delta,varargin)
%P = STBLPDF(X,ALPHA,BETA,GAM,DELTA) returns the pdf of the stable 
% distribtuion with characteristic exponent ALPHA, skewness BETA, scale
% parameter GAM, and location parameter DELTA, at the values in X.  We use 
% the parameterization of stable distribtuions used in [2] - The 
% characteristic function phi(t) of a S(ALPHA,BETA,GAM,DELTA)
% random variable has the form
%
% phi(t) = exp(-GAM^ALPHA |t|^ALPHA [1 - i BETA (tan(pi ALPHA/2) sign(t)]
%                  + i DELTA t )  if alpha ~= 1
%
% phi(t) = exp(-GAM |t| [ 1 + i BETA (2/pi) (sign(t)) log|t|] + i DELTA t
%                                 if alpha = 1
%
% The size of P is the size of X.  ALPHA,BETA,GAM and DELTA must be scalars
% 
%P = STBLPDF(X,ALPHA,BETA,GAM,DELTA,TOL) computes the pdf to within an
% absolute error of TOL.
%
% The algorithm works by computing the numerical integrals in Theorem
% 1 in [1] using MATLAB's QUADV function.  The integrands  
% are smooth non-negative functions, but for certain parameter values 
% can have sharp peaks which might be missed.  To avoid this, STBLEPDF
% locates the maximum of this integrand and breaks the integral into two
% pieces centered around this maximum (this is exactly the idea suggested
% in [1] ).  
%
% If abs(ALPHA - 1) < 1e-5,  ALPHA is rounded to 1.
%
%P = STBLPDF(...,'quick') skips the step of locating the peak in the 
% integrand, and thus is faster, but is less accurate deep into the tails
% of the pdf.  This option is useful for plotting.  In place of 'quick',
% STBLPDF also excepts a logical true or false (for quick or not quick)
%
% See also: STBLRND, STBLCDF, STBLINV, STBLFIT
%
% References:
%
% [1] J. P. Nolan (1997)
%     "Numerical Calculation of Stable Densities and Distribution
%     Functions"  Commun. Statist. - Stochastic Modles, 13(4), 759-774
%
% [2] G Samorodnitsky, MS Taqqu (1994)
%     "Stable non-Gaussian random processes: stochastic models with 
%      infinite variance"  CRC Press
%

if nargin < 5
    error('stblpdf:TooFewInputs','Requires at least five input arguments.'); 
end

% Check parameters
if alpha <= 0 || alpha > 2 || ~isscalar(alpha)
    error('stblpdf:BadInputs',' "alpha" must be a scalar which lies in the interval (0,2]');
end
if abs(beta) > 1 || ~isscalar(beta)
    error('stblpdf:BadInputs',' "beta" must be a scalar which lies in the interval [-1,1]');
end
if gam < 0 || ~isscalar(gam)
    error('stblpdf:BadInputs',' "gam" must be a non-negative scalar');
end
if ~isscalar(delta)
    error('stblpdf:BadInputs',' "delta" must be a scalar');
end

% Warn if alpha is very close to 1 or 0
if ( 1e-5 < abs(1 - alpha) && abs(1 - alpha) < .02) || alpha < .02 
    warning('stblpdf:ScaryAlpha',...
        'Difficult to approximate pdf for alpha close to 0 or 1')
end

% warnings will happen during call to QUADV, and it's okay
warning('off');

% Check and initialize additional inputs
quick = false;
tol = [];
for i=1:length(varargin)
    if strcmp(varargin{i},'quick')
        quick = true;
    elseif islogical(varargin{i})
        quick = varargin{end};
    elseif isscalar(varargin{i})
        tol = varargin{i};
    end
end

if isempty(tol)
    if quick 
        tol = 1e-8;
    else
        tol = 1e-12;
    end
end
        

%======== Compute pdf ==========%

% Check to see if you are in a simple case, if so be quick, if not do
% general algorithm
if alpha == 2                  % Gaussian distribution 
    x = (x - delta)/gam;                 % Standardize
    p = 1/sqrt(4*pi) * exp( -.25 * x.^2 ); % ~ N(0,2)
    p = p/gam; %rescale

elseif alpha==1 && beta == 0   % Cauchy distribution
    x = (x - delta)/gam;              % Standardize
    p = (1/pi) * 1./(1 + x.^2); 
    p = p/gam; %rescale

elseif alpha == .5 && abs(beta) == 1 % Levy distribution 
    x = (x - delta)/gam;              % Standardize
    p = zeros(size(x));
    if  beta ==1
        p( x <= 0 ) = 0;
        p( x > 0 ) = sqrt(1/(2*pi)) * exp(-.5./x(x>0)) ./...
                                                x(x>0).^1.5;
    else
        p(x >= 0) = 0;
        p(x < 0 ) = sqrt(1/(2*pi)) * exp(.5./x(x<0)  ) ./...
                                            ( -x(x<0) ).^1.5;
    end
    p = p/gam; %rescale
    
elseif abs(alpha - 1) > 1e-5          % Gen. Case, alpha ~= 1
    
    xold = x; % Save for later
    % Standardize in (M) parameterization ( See equation (2) in [1] ) 
    x = (x - delta)/gam - beta * tan(alpha*pi/2);  
    
    % Compute pdf
    p = zeros(size(x));
    zeta = - beta * tan(pi*alpha/2);  
    theta0 = (1/alpha) * atan(beta*tan(pi*alpha/2));
    A1 = alpha*theta0;
    A2 = cos(A1)^(1/(alpha-1));
    exp1 = alpha/(alpha-1);
    alpham1 = alpha - 1;
    c2 = alpha ./ (pi * abs(alpha - 1) * ( x(x>zeta) - zeta) ); 
    V = @(theta) A2 * ( cos(theta) ./ sin( alpha*(theta + theta0) ) ).^exp1.*...
        cos( A1 + alpham1*theta ) ./ cos(theta);
    
    
    % x > zeta, calculate integral using QUADV
    if any(x(:) > zeta)
        xshift = (x(x>zeta) - zeta) .^ exp1;
        
        if beta == -1 && alpha < 1
            p(x > zeta) = 0;
        elseif ~quick % Locate peak in integrand and split up integral        
            g = @(theta) xshift(:) .* V(theta) - 1;
            R = repmat([-theta0, pi/2 ],numel(xshift),1);
            if abs(beta) < 1
                theta2 = bisectionSolver(g,R,alpha);
            else
                theta2 = bisectionSolver(g,R,alpha,beta,xshift);
            end
            theta2 = reshape(theta2,size(xshift));
            % change variables so the two integrals go from 
            % 0 to 1/2 and 1/2 to 1.
            theta2shift1 = 2*(theta2 + theta0);
            theta2shift2 = 2*(pi/2 - theta2);
            g1 = @(theta)  xshift .* ...
                V(theta2shift1 * theta - theta0);
            g2 = @(theta)  xshift .* ...
                V(theta2shift2 * (theta - .5) + theta2);
            zexpz = @(z) max(0,z .* exp(-z)); % use max incase of NaN
           
            p(x > zeta) = c2 .* ...
                (theta2shift1 .* quadv(@(theta) zexpz( g1(theta) ),...
                                        0 , .5, tol) ...
               + theta2shift2 .* quadv(@(theta) zexpz( g2(theta) ),...
                                       .5 , 1, tol) );                       
                              
        else  % be quick - calculate integral without locating peak
              % Use a default tolerance of 1e-6
            g = @(theta) xshift * V(theta);
            zexpz = @(z) max(0,z .* exp(-z)); % use max incase of NaN
            p( x > zeta ) = c2 .* quadv(@(theta) zexpz( g(theta) ),...
                                        -theta0 , pi/2, tol );  
        end
        p(x > zeta) = p(x>zeta)/gam; %rescale
        
    end
    
    % x = zeta, this is easy
    if any( abs(x(:) - zeta) < 1e-8 )  
        p( abs(x - zeta) < 1e-8 ) = max(0,gamma(1 + 1/alpha)*...
            cos(theta0)/(pi*(1 + zeta^2)^(1/(2*alpha))));
        p( abs(x - zeta) < 1e-8 ) = p( abs(x - zeta) < 1e-8 )/gam; %rescale
        
    end
   
    % x < zeta, recall function with -xold, -beta, -delta 
    % This doesn't need to be rescaled.
    if any(x(:) < zeta)
        p( x < zeta ) = stblpdf( -xold( x<zeta ),alpha,-beta,...
                        gam , -delta , tol , quick); 
    end
        
else                    % Gen case, alpha = 1
    
    x = (x - (2/pi) * beta * gam * log(gam) - delta)/gam; % Standardize
    
    % Compute pdf
    piover2 = pi/2;
    twooverpi = 2/pi;
    oneoverb = 1/beta;
    theta0 = piover2;
    % Use logs to avoid overflow/underflow
    logV = @(theta) log(twooverpi * ((piover2 + beta *theta)./cos(theta))) + ...
                 ( oneoverb * (piover2 + beta *theta) .* tan(theta) );
    c2 = 1/(2*abs(beta));
    xterm = ( -pi*x/(2*beta));
    
    if ~quick  % Locate peak in integrand and split up integral
             % Use a default tolerance of 1e-12
        logg = @(theta) xterm(:) + logV(theta) ;
        R = repmat([-theta0, pi/2 ],numel(xterm),1);
        theta2 = bisectionSolver(logg,R,1-beta);     
        theta2 = reshape(theta2,size(xterm));
        % change variables so the two integrals go from 
        % 0 to 1/2 and 1/2 to 1.
        theta2shift1 = 2*(theta2 + theta0);
        theta2shift2 = 2*(pi/2 - theta2);
        logg1 = @(theta)  xterm + ...
            logV(theta2shift1 * theta - theta0);
        logg2 = @(theta)  xterm + ...
            logV(theta2shift2 * (theta - .5) + theta2);
        zexpz = @(z) max(0,exp(z) .* exp(-exp(z))); % use max incase of NaN

        p = c2 .* ...
            (theta2shift1 .* quadv(@(theta) zexpz( logg1(theta) ),...
                                    0 , .5, tol) ...
           + theta2shift2 .* quadv(@(theta) zexpz( logg2(theta) ),...
                                   .5 , 1, tol) );     
      
       
    else % be quick - calculate integral without locating peak
              % Use a default tolerance of 1e-6
        logg = @(theta) xterm + logV(theta);
        zexpz = @(z) max(0,exp(z) .* exp(-exp(z))); % use max incase of NaN
        p = c2 .* quadv(@(theta) zexpz( logg(theta) ),-theta0 , pi/2, tol );
            
    end
    
    p = p/gam; %rescale
    
end

p = real(p); % just in case a small imaginary piece crept in 
             % This might happen when (x - zeta) is really small   

end




function X = bisectionSolver(f,R,alpha,varargin)
% Solves equation g(theta) - 1 = 0 in STBLPDF using a vectorized bisection 
% method and a tolerance of 1e-5.  The solution to this
% equation is used to increase accuracy in the calculation of a numerical
% integral.   
%
% If alpha ~= 1 and 0 <= beta < 1, the equation always has a solution
%
% If alpha > 1 and beta <= 1, then g is monotone decreasing
% 
% If alpha < 1 and beta < 1, then g is monotone increasing
%
% If alpha = 1,  g is monotone increasing if beta > 0 and monotone 
% decreasing is beta < 0.  Input alpha = 1 - beta to get desired results.
%
%


if nargin < 2
    error('bisectionSolver:TooFewInputs','Requires at least two input arguments.'); 
end

noSolution = false(size(R,1));
% if ~isempty(varargin)
%     beta = varargin{1};
%     xshift = varargin{2};
%     if abs(beta) == 1
%         V0=(1/alpha)^(alpha/(alpha-1))*(1-alpha)*cos(alpha*pi/2)*xshift;
%         if alpha > 1
%             noSolution = V0 - 1 %>= 0;
%         elseif alpha < 1
%             noSolution = V0 - 1 %<= 0;
%         end
%     end 
% end
    
tol = 1e-6;
maxiter = 30;
    
[N M] = size(R);
if M ~= 2
    error('bisectionSolver:BadInput',...
        '"R" must have 2 columns');
end

a = R(:,1);
b = R(:,2);
X = (a+b)/2;

try
    val = f(X);
catch ME
    error('bisectionSolver:BadInput',...
        'Input function inconsistint with rectangle dimension')
end
  
if size(val,1) ~= N
    error('bisectionSolver:BadInput',...
        'Output of function must be a column vector with dimension of input');
end

% Main loop
val = inf;
iter = 0;

while( max(abs(val)) > tol && iter < maxiter )
    X = (a + b)/2;
    val = f(X);
    l = (val > 0);
    if alpha > 1
        l = 1-l;
    end
    a = a.*l + X.*(1-l);
    b = X.*l + b.*(1-l);
    iter = iter + 1;
end



if any(noSolution(:))
    X(noSolution) = (R(1,1) + R(1,2))/2;
end

end














