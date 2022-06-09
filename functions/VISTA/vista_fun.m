function [PEInd, FRInd, samp] = vista_fun(param)
% Author: Rizwan Ahmad (ahmad.46@osu.edu)

% Output
% ky(i) = PEInd(i), 'i' is the order of acquisition
% t(i) = FRInd(i) 
% samp = binary mask in ky-t domain

n     = param.n;     % Number of samples per frame
N     = param.PE;     % Number of phase encoding steps
FR    = param.FR;     % Number of frames
tf    = param.tf;    % Relative step-size in "time" dimension wrt to phase-encoding direction
w     = param.w;     % Scaling of time dimension
fs    = param.fs;    % To make time average fully sampled or not
uni   = param.uni;   % Iteration where sampling is reset to jittered uniform
R     = param.PE / param.n;     % Net acceleration factor
nIter = param.nIter; % Number of iterations
s     = log10(param.s);  % 1=uniform density, s>1 for variable density
sig   = param.sig;   % Std of the Gaussian envelope that defines variable density
beta  = param.beta;   % Exponent of force term
g     = param.g;     % At every g^th iteration, sampled are relocated to the nearest Cartesian grid
sz    = param.sz;    % Size of displayed samples
fl    = param.fl;    % At/beyond this iteration, start checking for fully-sampledness
ss    = param.ss;    % Gradient descent step-size
dsp   = param.dsp;   % Display the distribution at every dsp^th iteration
sd    = param.sd;    % Seed for random number generation
M     = param.M; % Total number of samples

%% Let's handle the special case of R = 1
if R == 1
    samp = ones(N, FR);
    return;
end


%% Use radom variable density sampling as initialization for VISTA
p1=(-floor(N/2):ceil(N/2)-1)';
t1=[];
Ind=0;
ti = zeros(n*FR,1);
ph = zeros(n*FR,1);
prob = (0.1 + s/(1-s+1e-10)*exp(-(p1).^2./(1*sig.^2)));

rng(sd);
tmpSd = round(1e6*rand(FR,1)); % Seeds for random numbers
for i = -floor(FR/2):ceil(FR/2)-1
    a = find(t1==i);
    n_tmp = n-numel(a);
    prob_tmp = prob;
    prob_tmp(a) = 0;
    p_tmp = randp(prob_tmp, tmpSd(i+floor(FR/2)+1), n_tmp, 1) - floor(N/2)-1;
    ti(Ind+1:Ind+n_tmp) = i;
    ph(Ind+1:Ind+n_tmp) = p_tmp;
    Ind = Ind+n_tmp;
end


% Displacement parameters
fprintf(['Computing VISTA, plese wait as it may take a while ...' '\n']);
tic,
stp = ones(1,nIter); % Gradient descent displacement
a = w*ones(1,nIter);  % Temporal axis scaling


f = 1+round(100*rand); % Figure index
dis_ext = zeros(M,1); % Extent of displacement
for i=1:nIter
    [ph,ti] = tile(ph(1:M), ti(1:M), param);
    for j=1:M
        
       
        % Distances -------------------------------------------------------
        dis= sqrt(abs(ph-ph(j)).^2 + abs(a(i)*(ti-ti(j))).^2);
        nanloc= dis==0;
        dis(nanloc)=inf;
        
        % Scaling ---------------------------------------------------------
        scl= 1 - s*exp(-(ph).^2./(2*sig.^2)); 
        scl = scl + (1-scl(1));
        dscl = 1/sig^2 * s.*ph(j).*exp(-(ph(j).^2)./(2*sig^2)); % Differentiation of scl wrt to "ph"

        % Force and resulting displacement --------------------------------
        fx = beta*(     (ph(j)-ph) .* (scl(j)*scl./dis.^(beta+2))) - dscl .* scl./dis.^beta;
        fy = beta*(a(i)^2*(ti(j)-ti) .* (scl(j)*scl./dis.^(beta+2)))*tf;
        ph(j) = (ph(j)+ max(min(stp(i)*sum(fx), R/4),-R/4));
        ti(j) = (ti(j)+ max(min(stp(i)*sum(fy), R/4),-R/4));
                
        % Ensure that the samples stay in bounds --------------------------
        if ph(j) < -floor(N/2)-1/2 
            ph(j) = ph(j) + N;
        elseif ph(j)> (ceil(N/2)-1/2)
            ph(j) = ph(j) - N;
        end

        if ti(j) < -floor(FR/2)-1/2 
            ti(j) = ti(j) + FR;
        elseif ti(j)> (ceil(FR/2)-1/2)
            ti(j) = ti(j) - FR;
        end
        
        % Displacing samples to nearest Cartesian location
        if rem(i,g)==0 || i==nIter
            ph(j)=round(ph(j));
            ti(j)=round(ti(j));
        end
        
        % Measuing the displacement
        if i==2, dis_ext(j) = abs(stp(i)*sum(fx)); end
    end
    
    
    % Normalizing the step-size to a reasonable value
    if i==3
        stp = ss*(1+R/4)/median(dis_ext)*stp;
    end
    
    % At uni-th iteration, reset to jittered uniform sampling
    ti = ti(1:M);
    ph = ph(1:M);
    if i == uni
        tmp = zeros(n, FR);
        for k = 1:FR
            tmp(:,k) = sort(ph((k-1)*n+1:k*n));
        end
        tmp = round(mean(tmp,2)); % Find average distances between adjacent phase-encoding samples
        ph = repmat(tmp, [FR,1]); % Variable density sampling with "average" distances
        rng(sd);
        rndTmp = rand(FR,1);
        for k=-floor(FR/2):ceil(FR/2)-1
            tmp = ti==k;
            ptmp = ph(tmp) + round(1/2*R^(1)*(rndTmp(k+floor(FR/2)+1)-0.5)); % Add jitter
            ptmp(ptmp>ceil(N/2)-1) = ptmp(ptmp>ceil(N/2)-1)-N; % Ensure the samples don't run out of the k-t grid
            ptmp(ptmp<-floor(N/2)) = ptmp(ptmp<-floor(N/2))+N;
            ph(tmp) = ptmp;
        end
        % Temporarily stretch the time axis to avoid bad local minima
        a(i+1:end) = a(i+1:end) .* (1 + exp(-((i+1:nIter)-(i+1))/ceil(nIter/60))); 
        %figure; plot(a); 
    end
    
    % Displace the overlapping points so that one location has only one sample
    if rem(i,g)==0 || i==nIter
        [ph, ti] = dispdup(ph(1:M), ti(1:M), param);
    end
    
    % Check/ensure time-average is fully sampled
    if (rem(i,g)==0 || i==nIter) && i>=fl 
        ph = ph(1:M);
        ti = ti(1:M);
        if fs==1 % Ensuring fully sampledness at average all
            [ph, ti] = fillK(ph, ti, ph, ti, param);
        elseif fs>1 % Ensuring fully sampledness for "fs" frames
            for m = 1:floor(FR/fs)
                tmp = (m-1)*n*fs + 1 : m*n*fs;
                [ph, ti] = fillK(ph(tmp), ti(tmp), ph, ti, param);
            end
        end
    end    
    
    if ((i==1 || rem(i,dsp)==0 || i==nIter) && dsp>0) % When to diplay the distribution
        figure(f);
        plot(ti(1:M),ph(1:M),'s','markersize',sz,'color','white','markerfacecolor','white'); 
        axis('image');axis([-floor(FR/2), ceil(FR/2)-1, -floor(N/2), ceil(N/2)-1]);
        title(['Iter: ' num2str(i) ', Number of samples: ' num2str(M)]);
        xlabel('t');
        ylabel('k_y');
        set(gca,'color','k')
        hold off;
    end
end
[ph, ti] = dispdup(ph(1:M), ti(1:M), param);

% create binary mask 'samp'
PEInd = round(N*(ti + floor(FR/2)) + (ph + floor(N/2)+1)); % shift center to N/2
samp=zeros(N, FR);
samp(PEInd)=1;

% flip the ky ordering of alternate frames to reduce jumps
FRind = size(PEInd);
[row, col] = ind2sub([N,FR], PEInd);
PEInd = row;
FRInd = col;
% PEind = zeros(M,1);
for i=1:FR
    if mod(i, 2) == 0
        % Sort in descending order (even frames)
        PEInd((i-1)*n+1 : i*n) = sort(PEInd((i-1)*n+1 : i*n), 'descend');
    else
        PEInd((i-1)*n+1 : i*n) = sort(PEInd((i-1)*n+1 : i*n), 'ascend');
    end
    FRind((i-1)*n+1 : i*n) = i;
end

fprintf('VISTA computed in %4.2f s\n', toc);

% Plotting
if dsp > 0
    E = 1; % Vista only supports one encoding
    tiFont = 20; % title font
    axFont = 14; % axis font
    laFont = 18; % label font
    figure;
    subNum = E + double((E/2) >= 1); % Number of subplots
    tiledlayout(1,subNum,'TileSpacing','compact', 'Padding', 'compact')
    for e = 1 : subNum
        nexttile;
%         subplot(1,subNum,e);
        if e == 1
            imagesc(samp(:,:,e),[0,E]); axis('image'); colormap(hot);
            set(gca, 'FontSize', axFont, 'FontName','times');
            xlabel('$t$', 'FontSize', laFont,'Interpreter','latex'); 
            ylabel('$k_y$', 'FontSize', laFont,'Interpreter','latex');
            if E==1 % do nothing
            elseif E>1, title(['Encoding ' num2str(e)], 'FontSize', tiFont,'Interpreter','latex'); end
        elseif e < (E+1)
            imagesc(samp(:,:,e),[0,E]); axis('image'); colormap(hot);
            set(gca, 'FontSize', axFont, 'FontName','times','ytick',[]);
            xlabel('$t$', 'FontSize', laFont,'Interpreter','latex');    
            title(['Encoding ' num2str(e)], 'FontSize', tiFont,'Interpreter','latex');
        else
            imagesc(max(samp,[],3)); axis('image'); colormap(hot);
            set(gca, 'FontSize', axFont, 'FontName','times','ytick',[]);
            xlabel('$t$', 'FontSize', laFont,'Interpreter','latex');  
            title('Superimposed', 'FontSize', tiFont,'Interpreter','latex');
        end
    end
    set(gcf,'color','w','units','points','position',[10,10,450,350]); %export_fig('vista',gcf,'-m4','-png');
    
    figure;
    len = min([length(PEInd), 120]);
    plot(1:len, PEInd(1:len,1), '-'); hold on;
    plot(1:len, PEInd(1:len,1), '.','MarkerSize', 12);
    set(gca, 'FontSize', axFont, 'FontName','times');
    xlabel('Acquisition Order', 'FontSize', laFont,'Interpreter','latex');
    ylabel('$k_y$', 'FontSize', laFont,'Interpreter','latex');
    title('PE Index vs. Acquisition Order', 'FontSize', tiFont,'Interpreter','latex');
    set(gcf,'color','w','units','points','position',[10,10,600,350]); %export_fig('vista-index',gcf,'-m4','-png');
end

samp = logical(samp); % convert to logical before returning


function [po,to] = tile(ph,ti,param)
% Replicate the sampling pattern in each direction. Probablity, this is 
% not an efficient way to impose periodic boundary condition because it 
% makes the problem size grow by 9 fold.

N = param.PE; 
FR = param.FR; 
po=[ph; ph-N; ph; ph+N;     ph-N; ph+N;    ph-N; ph; ph+N];
to=[ti; ti-FR; ti-FR; ti-FR;   ti; ti;        ti+FR; ti+FR; ti+FR];

