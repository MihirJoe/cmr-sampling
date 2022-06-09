function[PEInd, FRInd, samp] = cava_fun(param)
% Author: Rizwan Ahmad (ahmad.46@osu.edu)

% Output
% ky(i,e) = PEInd(i,e), where 'e' is encoding and 'i' is the order of acquisition
% t(i) = FRInd(i)
% samp = binary mask in ky-t-encoding domain

n   = param.n;    % Number of phase encoding (PE) lines per frame. This can be changed post-acquisition
M   = param.M; % Total number of samples
N   = param.PE;   % Size of of PE grid
E   = param.E;    % Number of encoding, E=1 for cine, E=2 for flow
tau = param.tau;
s   = param.s; 
a   = param.alph;
dsp = param.dsp;


R  = N/n; % Initial guess of acceleration
% s  = max(1, s*(R).^(1/3));  % This is to adapt the value of's' based on acceleration.

gr = (1+sqrt(5))/2; % golden ratio F_{PE+1}/F_{PE}
gr = 1/(gr+tau-1); %(1-1/gr); % golden angle 
Ns = ceil(N * 1/s); % Size of shrunk PE grid

%% Size of the smaller pseudo-grid which after stretching gives the true grid size
k = (N/2-Ns/2)/((Ns/2)^a); %(1/2*(N-Ns)/max(tmp)); % location specific displacement


%% Let's populate the grid;
samp  = zeros(N, ceil(M/n), E); % sampling on k-t grid
PEInd = zeros(M,E); % Final PE index used for MRI
FRInd = zeros(M,1); % Frame index (can be changed post acquisition)

ind = zeros(M,E); % "hidden" index on a uniform grid

for e=1:(E+1)
    %t = tiledlayout(1, E, 'TileSpacing', 'Compact');
    kk=e;
    for i=1:M
        if e < E + 1
            if i==1,     ind(i,e) = rem(floor(Ns/2) + 1 + (e-1)*sqrt(11)*gr*Ns/E -1, Ns) + 1;
            elseif i>1,  ind(i,e) = rem((ind(i-1,e) + Ns*gr)-1, Ns) + 1;
            end
            ind(i,e) = ind(i,e) - Ns.*(ind(i,e)>=(Ns+0.5));
            
            if rem(N,2)==0 % if even, shift by 1/2 pixel
                indC = ind(i,e) - k*sign((Ns/2+1/2)-ind(i,e))*(abs((Ns/2+1/2)-ind(i,e))).^a + (N-Ns)/2 + 1/2;%(ctrn - ctrnS);
                indC = indC - N.*(indC>=(N+0.5));
            elseif rem(N,2)==1 % if odd don't shift
                indC = ind(i,e) - k*sign((Ns/2+1/2)-ind(i,e))*(abs((Ns/2+1/2)-ind(i,e))).^a + (N-Ns)/2;%(ctrn - ctrnS);
            end
            %     kk = (-1)^i;
            PEInd(i,e) = round(indC);
            samp(PEInd(i,e), ceil(i/n),e) = samp(PEInd(i,e), ceil(i/n),e)+ kk;
        end
        FRInd(i) = ceil(i/n);
    end
end

if dsp == 1
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
            title(['Encoding ' num2str(e)], 'FontSize', tiFont,'Interpreter','latex');
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
    set(gcf,'color','w','units','points','position',[10,10,450,350]); %export_fig('cava',gcf,'-m4','-png');
    
    figure;
    len = min([length(PEInd), 120]);
    plot(1:len, PEInd(1:len,1), '-'); hold on;
    plot(1:len, PEInd(1:len,1), '.','MarkerSize', 12);
    set(gca, 'FontSize', axFont, 'FontName','times');
    xlabel('Acquisition Order', 'FontSize', laFont,'Interpreter','latex');
    ylabel('$k_y$', 'FontSize', laFont,'Interpreter','latex');
    title('PE Index vs. Acquisition Order', 'FontSize', tiFont,'Interpreter','latex');
    set(gcf,'color','w','units','points','position',[10,10,600,350]); %export_fig('cava-index',gcf,'-m4','-png');
end

samp = logical(samp); % convert to logical before returning
