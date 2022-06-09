function [PEInd] = opra_fun(param)

% Output
% ky(i) = PEInd(i,1); kz(i) = PEInd(i,2) where 'i' is the order of acquisition

N   = param.PE; % matrix size
n   = param.n; % Nominal number of samples in each frame (display only)
FR  = param.FR; % Nominal number of frames (display only)
M   = param.M; % Total number of samples
L   = param.L; % size of leaflet
s   = param.s/2; % controls sampling density
ar  = param.ar+1; % ar>=0; controls the aspect ratio/orientation of the high density ellipsoid region in the middle 
ga  = (sqrt(5)-1)/2 * pi; % golden angle
gs  = param.gs;% second irrational number
cg  = param.cg; %size of the gap in the middle
phi = param.phi; % Angular jump from the end of one leaflet to the start of the next
dsp = param.dsp; % plot figures or not
p  = 0; % A larger value avoids a large number of samples at the center. The sampling density is not sensitity to this parameter


 
th0 = 0;
N_0 = N; % original N
N = N/2; % radii of the sampling grid
x = nan(ceil(M/L)*L, 1); % x index
y = nan(ceil(M/L)*L, 1); % y index
% avgSamp = zeros(2*N); % time average sampling 
PEInd = zeros(ceil(M/L)*L,2); % sampling indices
% figure; 
for i=1:ceil(M/L)
    tSamp = zeros(2*N); % temporary variable
    eAng = ga - phi; % angle of the elbow
    
    % first half of "L"
    th1 = rem(th0 + (i-1)*ga, 2*pi);
    th1 = atan2(sin(th1)*(N(2)/N(1))^ar, cos(th1));
    R1 = (N(1)*N(2))/sqrt((N(1)*sin(th1))^2 + (N(2)*cos(th1))^2); % radius for first half of "L"
    del1 = 2*R1/L; % Average distance between two samples
    r1 = (R1:-del1:R1/(L/2)) - rem((i-1)*gs,1)*del1; 
    r1 = (r1+p).^s;
    r1 = (r1/((R1+p)^(s-1))) * (R1-1-cg-p)/R1 + cg;
    
    % second half of "L"
    th2 = rem(th0 + (i-1)*ga + eAng, 2*pi);
    th2 = atan2(sin(th2)*(N(2)/N(1))^ar, cos(th2));
    R2 = (N(1)*N(2))/sqrt((N(1)*sin(th2))^2 + (N(2)*cos(th2))^2); % radius for second half of "L"
    del2 = 2*R2/L;
    r2 = (0:del2:R2-R2/(L/2)) + rem((i-1)*gs,1)*del2;
    r2 = (r2+p).^s;
    r2 = (r2/((R2+p)^(s-1))) * (R2-1-cg-p)/R2 + cg;

    % polar to Cartesian
    x1 = r1*cos(th1)+N(1)+1-1e-6;
    y1 = r1*sin(th1)+N(2)+1-1e-6;
    x2 = r2*cos(th2)+N(1)+1-1e-6;
    y2 = r2*sin(th2)+N(2)+1-1e-6;
    x(1+(2*i-2)*L/2:(2*i-1)*L/2) = x1;
    x(1+(2*i-1)*L/2:(2*i-0)*L/2) = x2;
    y(1+(2*i-2)*L/2:(2*i-1)*L/2) = y1;
    y(1+(2*i-1)*L/2:(2*i-0)*L/2) = y2;
end
x = round(x); x(x > N_0(1)) = N_0(1); x(x<1) = 1;
y = round(y); y(y > N_0(2)) = N_0(2); y(y<1) = 1;

PEInd(:,1) = x; % sampling index
PEInd(:,2) = y;

% Plotting
if dsp == 1
    tiFont = 20; % title font
    axFont = 14; % axis font
    laFont = 18; % label font
    figure;
    rang = [1, n]; % the number of samples plot in each frame
    col = 5;
    tiledlayout(1,col,'TileSpacing','compact', 'Padding', 'compact')
    for i=1:col % position in subplot space
        nexttile;
        if i < 4
            %Plot samples
            tSamp = zeros(N_0);
            for j=rang(1):rang(2)
              tSamp(PEInd(j,1), PEInd(j,2)) = tSamp(PEInd(j,1), PEInd(j,2)) + 1;
            end
            imagesc(sqrt(tSamp)); axis('image'); colormap(gray);
            set(gca, 'FontSize', axFont, 'FontName','times');
            xlabel('$k_z$', 'FontSize', laFont,'Interpreter','latex'); 
            if i == 1, ylabel('$k_y$', 'FontSize', laFont,'Interpreter','latex');
            else i>1 && i <4; set(gca, 'FontSize', axFont, 'FontName','times','ytick',[]); end
            title([num2str(rang(1)) '-' num2str(rang(2)) ' Samples'], 'FontSize', tiFont,'Interpreter','latex');
            
            rang(1) = rang(2) + 1;
            rang(2) = rang(2) + n;
        
        elseif i == 4
            %Plot accumulation
            tSamp = zeros(N_0);
            for j=1:M
              tSamp(PEInd(j,1), PEInd(j,2)) = tSamp(PEInd(j,1), PEInd(j,2)) + 1;
            end
            imagesc(tSamp.^0.25); axis('image'); colormap(gray);
            set(gca, 'FontSize', axFont, 'FontName','times','ytick',[]);
            xlabel('$k_z$', 'FontSize', laFont,'Interpreter','latex'); 
            title([num2str(1) '-' num2str(M) ' Samples'], 'FontSize', tiFont,'Interpreter','latex');
          
        elseif i == 5
            %Plot binary mask
            imagesc(logical(tSamp)); axis('image'); colormap(gray);
            set(gca, 'FontSize', axFont, 'FontName','times','ytick',[]);
            xlabel('$k_z$', 'FontSize', laFont,'Interpreter','latex'); 
            title(['Binary Mask'], 'FontSize', tiFont,'Interpreter','latex');
        end
    end
    set(gcf,'color','w','units','points','position',[10,10,900,300]); %export_fig('opra',gcf,'-m4','-png');
end