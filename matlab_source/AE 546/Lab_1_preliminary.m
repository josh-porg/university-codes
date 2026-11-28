clear; close; clc;

%bring in functions from ae 445
addpath(genpath('AE 445'));


% Airfoil Geometry
airfoil = loadAirfoil();
[DxLower_c, DxUpper_c, DyLower_c, DyUpper_c, xLower_c, xUpper_c, yLower_c, yUpper_c] = possitionDeltas(airfoil)

% Experimental Data
opts = detectImportOptions("Lab1_Data.xlsx");
dataMatrix = readmatrix("Lab1_Data.xlsx", "Range", "A1:G51"); % Data in A2 trhough G51
% first row is angle of attack second row is static head thid row is total
% head subsiquent rows are local head measuremnts
AoAs = dataMatrix(1,:);
staticHeads = dataMatrix(2,:);
totalHeads = dataMatrix(3,:);
localHeads = dataMatrix(4:end,:); % I belice all is well with data imports

% initialize constants
rhoAir = 1.225; % kg/m^3
rhoLiquid = 826; % kg/m^3
c = .195; % m


%variabel for data to be saved in
countOfAoAs = length(AoAs);
countOPoints = length(localHeads);
Res = zeros(0, countOfAoAs);
Cls = zeros(0, countOfAoAs);
Cds = zeros(0, countOfAoAs);
Cmc_4s = zeros(0, countOfAoAs);
CpUppers  = zeros(countOPoints / 2, countOfAoAs);
CpLowers  = zeros(countOPoints / 2, countOfAoAs);


%TODO: repeat for each AOA
i=1 % make i for loop itteration var
for i = 1:countOfAoAs
    
    alpha = AoAs(i)
    hs = staticHeads(i);
    ht = totalHeads(i);
    hlocal = localHeads(:,i);
    %[CpUpper, CpLower] = getCps(hs, ht, hlocal)
    [CpUppers(:,i), CpLowers(:,i)] = getCps(hs, ht, hlocal)
    
    vInf = velocityInf(ht,hs)
    
    % clalculate Rn for an AoA
    %Re = reynoldsNumber(vInf)
    Res(i) = reynoldsNumber(vInf)
    
    % calculate at one AoA; local coeffs: cp cl cd cm,c/4
    %[Cl, Cd, CmC_4] = getCoeffs(alpha, CpUpper, CpLower, DxLower_c, DxUpper_c, DyLower_c, DyUpper_c,xLower_c, xUpper_c, yLower_c, yUpper_c)
    % [Cls(i), Cds(i), Cmc_4s(i)] = getCoeffs(alpha, CpUpper, CpLower, DxLower_c, DxUpper_c, DyLower_c, DyUpper_c,xLower_c, xUpper_c, yLower_c, yUpper_c)
    [Cls(i), Cds(i), Cmc_4s(i)] = getCoeffs(alpha, CpUppers(i), CpLowers(i), DxLower_c, DxUpper_c, DyLower_c, DyUpper_c,xLower_c, xUpper_c, yLower_c, yUpper_c)
    %TODO: I belive there to be an error in get coeefs becuase cl of 1 seems
    %too high
    
    % a bunch of other stuff
end


% Make charts for forty-eight distributions of pressure coefficients corresponding to all AoA s. 
% Indicate the Reynolds number of each AoA. 
% Having graphical illustrations would be better, but not mandatory.

figure;
plotScalarVsAoA([AoAs; Cls; Cds; Cmc_4s])
figure;
plotScalarVsCl(Cls, Cds, Cmc_4s)
figure;
plotReVsAoa(Res,AoAs)

%probably put this inside the loop
for i = 1:countOfAoAs
    figure;
    %dont forget to put Re in these plots
    plotCpDistribution(airfoil,CpUppers(:,i),CpLowers(:,i),AoAs(i))
    figure;
    plotCpDistributionWithAirfoil(airfoil,CpUppers(:,i),CpLowers(:,i),AoAs(i))
end

%Make figures for total coefficients of lift, drag and moment on the quarter chord corresponding to different AoA s, just like Fig.1&2.
% These results should be compared  with theoretical results, which can either be calculated by the classical thin airfoil theory or extracted from Fig.1 & 2


function Re = reynoldsNumber(v)
    mu = 1.78E-5; % kg/(m*s)
    rho = 1.225; % kg/m^3
    c = .195; % m
    Re = rho * v * c / mu;
end

function vInf = velocityInf(ht,hs)
    rhoAir = 1.225; % kg/m^3
    rhoLiquid = 826; % kg/m^3
    g = 9.81; % m/s^2
    vInf = sqrt(2 * rhoLiquid * g * abs(ht-hs) / rhoAir); % m/s
end

function airfoil = loadAirfoil()
    iaf.designation='23012';
    iaf.n=24;
    iaf.HalfCosineSpacing=1;
    iaf.wantFile=1;
    iaf.datFilePath='./'; % Current folder
    iaf.is_finiteTE=0;
    airfoil = naca5gen(iaf);
end

function [DxLower_c, DxUpper_c, DyLower_c, DyUpper_c, xLower_c, xUpper_c, yLower_c, yUpper_c] = possitionDeltas(airfoil)
    % I belive this need updating and position os sensors is average of the
    % beteen points given
    xL = airfoil.xL;
    xU = airfoil.xU;
    yL = airfoil.zL;
    yU = airfoil.zU;
    %careful with bounds runfrom begining+1 to end-1
    DxLower_c = zeros(length(xL),0);
    DxUpper_c = zeros(length(xL),0);
    DyLower_c = zeros(length(xL),0);
    DyUpper_c = zeros(length(xL),0);
    
    xLower_c = zeros(length(xL),0);
    xUpper_c = zeros(length(xL),0);
    yLower_c = zeros(length(xL),0);
    yUpper_c = zeros(length(xL),0);
    
    for i = 2:length(xL)-1
        DxLower_c(i) = abs(1/2 * (xL(i+1) - xL(i-1)));
        DxUpper_c(i) = abs(1/2 * (xU(i+1) - xU(i-1)));
        DyLower_c(i) = 1/2 * (yL(i+1) - yL(i-1));
        DyUpper_c(i) = 1/2 * (yU(i+1) - yU(i-1));
        
%         xLower_c(i) = xL(i); %TODO: should this be + the average?
%         xUpper_c(i) = xU(i);
%         yLower_c(i) = yL(i);
%         yUpper_c(i) = yU(i);
    end
        %testing this out
        xLower_c = ((xL(1:end-1) + xL(2:end)) / 2)'; %TODO: should this be + the average?
        xUpper_c = ((xU(1:end-1) + xU(2:end)) / 2)'; 
        yLower_c = ((yL(1:end-1) + yL(2:end)) / 2)'; 
        yUpper_c = ((yU(1:end-1) + yU(2:end)) / 2)'; 
        
        %just a guess
        %DxLower_c = (xLower_c - xL(1:end-1)') * 2;
    
    
end

function cp = getCp(hs, ht, h)
    cp = (h-hs)./(ht-hs);
end

function [CpUpper, CpLower] = getCps(hs, ht, hlocal)
    hUpper = hlocal(1:end/2);
    hLower = hlocal(end/2+1:end);
    CpUpper = getCp(hs, ht, hUpper);
    CpLower = getCp(hs, ht, hLower);
end


function [Cn, Ca, CmLe] = coefficientsOfForce(CpUpper, CpLower, DxLower_c, DxUpper_c, DyLower_c, DyUpper_c, xLower_c, xUpper_c, yLower_c, yUpper_c)
    Cn = sum(CpLower) * sum(DxLower_c) - sum(CpUpper) * sum(DxUpper_c) ;
    Ca = sum(CpLower) * sum(DyLower_c) - sum(CpUpper) * sum(DyUpper_c);
    CmLe = sum(CpUpper) * sum(DxUpper_c) * sum(xUpper_c) - sum(CpLower) * sum(DxLower_c) * sum(xLower_c) + sum(CpUpper) * sum(DyUpper_c) * sum(yUpper_c) - sum(CpLower) * sum(DyLower_c) * sum(yLower_c);
end

function [Cl, Cd, CmC_4] = coefficientsOfAerodynamicForces(alpha, CoefForce)
    Cn = CoefForce(1);
    Ca = CoefForce(2);
    CmLe = CoefForce(3);
    Cl = Cn * cosd(alpha) - Ca * sind(alpha);
    Cd = Cn * sind(alpha) + Ca * cosd(alpha);
    CmC_4 = CmLe + Cl / 4;
end

function [Cl, Cd, CmC_4] = getCoeffs(alpha, CpUpper, CpLower, DxLower_c, DxUpper_c, DyLower_c, DyUpper_c, xLower_c, xUpper_c, yLower_c, yUpper_c)
    [Cn, Ca, CmLe] = coefficientsOfForce(CpUpper, CpLower, DxLower_c, DxUpper_c, DyLower_c, DyUpper_c, xLower_c, xUpper_c, yLower_c, yUpper_c);
    CoefForce = [Cn, Ca, CmLe];
    [Cl, Cd, CmC_4] = coefficientsOfAerodynamicForces(alpha, CoefForce);
end

function plotAirfoil(airfoil,alpha)
    % TODO: consider implementing alpha later
    % for alpha we want to move airfoil to center of rotation, rotate and
    % then translate it back
    % TODO: mayby implement alpha by passing a transformation matrix
    x = airfoil.x; z = airfoil.z;
    xU = airfoil.xU; zU = airfoil.zU;
    xL = airfoil.xL; zL = airfoil.zL;
    xC = airfoil.xC; zC = airfoil.zC;
    name = airfoil.name;
    
    hold on;
    airfoilPlot(1) = plot(x,z,'b-');
    set(airfoilPlot(1),'DisplayName','Airfoil Surface');
%     plot(xU,zU,'bo-')
%     plot(xL,zL,'ro-')

    airfoilPlot(2) = plot(xC,zC,'r--');
    set(airfoilPlot(2),'DisplayName','Mean Camber Line');

    title(name)
end

function plotCpDistribution(airfoil,CpUpper,CpLower,alpha)
    x = airfoil.x; z = airfoil.z;
    xU = airfoil.xU; zU = airfoil.zU;
    xL = airfoil.xL; zL = airfoil.zL;
    xC = airfoil.xC; zC = airfoil.zC;
    
    %consider plotting lift at that cord value as well
    % calculate pressure sensor locations
    xUs = (xU(1:end-1) + xU(2:end)) / 2;
    xLs = (xL(1:end-1) + xL(2:end)) / 2;
    
    plot1(1) = plot(xUs,CpUpper,'yo-'); hold on;
    plot1(2) = plot(xLs,CpLower,'co-'); hold on;
    
    set(plot1(1),'DisplayName','C_pUpper','Color',[1 1 0]); set(plot1(2),'DisplayName','C_pLower','Color',[0 1 1]);
    set(gca,'YDir','reverse'); %reversed because it makes sense as up and down pressures
    
    title(sprintf("C_p Distribution @ alpha = %g", alpha),"Color",'w'); xlabel('x / c ','Color','w'); ylabel('Coefficient of Pressure C_p','Color','w');
    
    set(gca, 'color','k', 'xcolor','w', 'ycolor','w'); set(gcf, 'color','k');
    grid(gca,"minor"); box(gca,'on'); hold(gca,'off');
    
    % Create legend
    legend1 = legend(gca,'show'); set(legend1,'TextColor',[1 1 1],'Orientation','vertical','Location','northeast','EdgeColor','none');
    
end

function plotCpDistributionWithAirfoil(airfoil,CpUpper,CpLower,alpha)
    x = airfoil.x; z = airfoil.z;
    xU = airfoil.xU; zU = airfoil.zU;
    xL = airfoil.xL; zL = airfoil.zL;
    xC = airfoil.xC; zC = airfoil.zC;
    
    %consider plotting lift at that cord value as well
    % calculate pressure sensor locations
    xUs = (xU(1:end-1) + xU(2:end)) / 2;
    xLs = (xL(1:end-1) + xL(2:end)) / 2;
    
    zUs = (zU(1:end-1) + zU(2:end)) / 2;
    zLs = (zL(1:end-1) + zL(2:end)) / 2;
    
%     plot1(1) = plot(xUs,-CpUpper + zUs,'yo-'); hold on;
%     plot1(2) = plot(xLs,-CpLower + zUs,'co-'); hold on;
    plot1(1) = plot(xUs,-CpUpper,'yo-'); hold on;
    plot1(2) = plot(xLs,-CpLower,'co-'); hold on;
    
    set(plot1(1),'DisplayName','C_pUpper','Color',[1 1 0]); set(plot1(2),'DisplayName','C_pLower','Color',[0 1 1]);
    %set(gca,'YDir','reverse');
    
    title(sprintf("C_p Distribution With airfoil @ alpha = %g", alpha),"Color",'w'); xlabel('x / c ','Color','w'); ylabel('Coefficient of Pressure C_p','Color','w');
    
    set(gca, 'color','k', 'xcolor','w', 'ycolor','w'); set(gcf, 'color','k');
    grid(gca,"minor"); box(gca,'on'); hold(gca,'off');
    
    % Create legend
    legend1 = legend(gca,'show'); set(legend1,'TextColor',[1 1 1],'Orientation','vertical','Location','northeast','EdgeColor','none');
    
    hold on
    plotAirfoil(airfoil,alpha)
    axis equal
    
    camroll(-alpha)
    hold on
    axis off
    
end

function plotReVsAoa(Re,AoA)
    plot1 = plot(AoA,Re,'wo-'); hold on;
    title("Reynolds Number V.S. Angle of Attack","Color",'w'); xlabel('Angle Of Attack, alpha, (deg)','Color','w'); ylabel('Reynolds Number','Color','w');
    set(plot1,'DisplayName','Reynolds Number','Color',[1 1 0]);
    set(gca, 'color','k', 'xcolor','w', 'ycolor','w'); set(gcf, 'color','k');
    grid(gca,"minor"); box(gca,'on'); hold(gca,'off');
end

function plotScalarVsAoA(scalars)
    AoA = scalars(1,:);
    %Re = scalars(2,:); %Re is aprox const
    Cl = scalars(2,:);
    Cd = scalars(3,:);
    Cmc_4 = scalars(4,:);
    
    plot1(1) = plot(AoA,Cl,'yo-'); hold on;
    plot1(2) = plot(AoA,Cd,'co-'); hold on;
    plot1(3) = plot(AoA,Cmc_4,'magentao-'); hold on;
    
    set(plot1(1),'DisplayName','C_l','Color',[1 1 0]); set(plot1(2),'DisplayName','C_d','Color',[0 1 1]); set(plot1(3),'DisplayName','C_m_,_c_/_4','Color',[1 0 1]);
    
    title("Coefficients V.S. Angle of Attack","Color",'w'); xlabel('Section Angle Of Attack, alpha, (deg)','Color','w'); ylabel('C_l, C_d, C_m_,_c_/_4','Color','w');
    
    set(gca, 'color','k', 'xcolor','w', 'ycolor','w'); set(gcf, 'color','k');
    grid(gca,"minor"); box(gca,'on'); hold(gca,'off');
    
    % Create legend
    legend1 = legend(gca,'show'); set(legend1,'TextColor',[1 1 1],'Orientation','vertical','Location','northeast','EdgeColor','none');
end

function plotScalarVsCl(Cl,Cd,Cmc_4)
    %Re = scalars(2,:); %Re is aprox const
%     Cl = scalars(1,:);
%     Cd = scalars(2,:);
%     Cmc_4 = scalars(3,:);
    
    plot1(1) = plot(Cl,Cd,'co-'); hold on;
    plot1(2) = plot(Cl,Cmc_4,'magentao-'); hold on;
    
    set(plot1(1),'DisplayName','C_d','Color',[0 1 1]); 
    set(plot1(2),'DisplayName','C_m_,_c_/_4','Color',[1 0 1]);
    
    title("C_d, C_m_,_c_/_4 V.S. C_l","Color",'w'); xlabel('Coefficient of Lift, C_l','Color','w'); ylabel('C_d, C_m_,_c_/_4','Color','w');
    
    set(gca, 'color','k', 'xcolor','w', 'ycolor','w'); set(gcf, 'color','k');
    grid(gca,"minor"); box(gca,'on'); hold(gca,'off');
    
    % Create legend
    legend1 = legend(gca,'show'); set(legend1,'TextColor',[1 1 1],'Orientation','vertical','Location','northeast','EdgeColor','none');
end
