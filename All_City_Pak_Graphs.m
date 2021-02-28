clear all 
clc
[Data] = xlsread("Data\Results\CSV\Energy Data.csv");
La = Data(:,1);
Lo = Data(:,2);
OT = Data(:,3);
FTR = Data(:,5);
FE = Data(:,6);
MTR = Data(:,7);
ME = Data(:,8);
Er = Data(:,9);
Hr = Data(:,10);
Er2 = Data(:,11);

%worldmap('pakistan')
%Pakistan=shaperead("E:\Books and pdfs\ISP\ShapeFile\PAK_adm1.shp")
%Pakistan1=geoshow(gca,[Pakistan.Y],[Pakistan.X],'linewidth',0.25);



% Optimal Tilt Radiation Map
figure(1);
latV = linspace(min(Lo), max(Lo), 100);
lonV = linspace(min(La), max(La), 100);
[latG, lonG] = meshgrid(latV, lonV);
valG = griddata(Lo, La, FTR, latG, lonG);
contourf(latG, lonG, valG, 100, 'LineColor', 'none')
xlabel('Longitude (°)')
ylabel('Latitude (°)')
colormap(jet)
a=colorbar;
ylabel(a,'Optimal Tilt Radiations (KW-hr / m^2 / Month)')
hold on;
bh = borders('pakistan','k-','LineWidth',2);
title('Optimal Tilt Irradiation')



%Optimal Energy
figure(2);
latV = linspace(min(Lo), max(Lo), 100);
lonV = linspace(min(La), max(La), 100);
[latG, lonG] = meshgrid(latV, lonV);
valG = griddata(Lo, La, FE, latG, lonG);
contourf(latG, lonG, valG, 100, 'LineColor', 'none')
xlabel('Longitude (°)')
ylabel('Latitude (°)')
colormap(jet)
a=colorbar;
ylabel(a,'Optimal Energy (KW-hr / m^2 / Month)')
hold on;
bh = borders('pakistan','k-','LineWidth',2);
title('Optimal Energy')



%Percentage Diff
figure(3);
latV = linspace(min(Lo), max(Lo), 100);
lonV = linspace(min(La), max(La), 100);
[latG, lonG] = meshgrid(latV, lonV);
valG = griddata(Lo, La, Er, latG, lonG);
contourf(latG, lonG, valG, 100, 'LineColor', 'none')
xlabel('Longitude (°)')
ylabel('Latitude (°)')
colormap(jet)
a=colorbar;
ylabel(a,'Percentage Difference %')
hold on;
bh = borders('pakistan','k-','LineWidth',2);
title('Monthly vs Fixed Energy Difference')


%Optimal Tilt
figure(4);
latV = linspace(min(Lo), max(Lo), 100);
lonV = linspace(min(La), max(La), 100);
[latG, lonG] = meshgrid(latV, lonV);
valG = griddata(Lo, La, OT, latG, lonG);
contourf(latG, lonG, valG, 100, 'LineColor', 'none')
xlabel('Longitude (°)')
ylabel('Latitude (°)')
colormap(jet)
a=colorbar;
ylabel(a,'Optimal Tilt (°)')
hold on;
bh = borders('pakistan','k-','LineWidth',2);
title('Optimal Tilt (°)')


% Horizontal Radiation Map
figure(5);
latV = linspace(min(Lo), max(Lo), 100);
lonV = linspace(min(La), max(La), 100);
[latG, lonG] = meshgrid(latV, lonV);
valG = griddata(Lo, La, Hr, latG, lonG);
contourf(latG, lonG, valG, 100, 'LineColor', 'none')
xlabel('Longitude (°)')
ylabel('Latitude (°)')
colormap(jet)
a=colorbar;
ylabel(a,'Global Horizontal Radiations (KW-hr / m^2 / Month)')
hold on;
bh = borders('pakistan','k-','LineWidth',2);
title('Global Horizontal Irradiation')


%Percentage Error
figure(6);
latV = linspace(min(Lo), max(Lo), 100);
lonV = linspace(min(La), max(La), 100);
[latG, lonG] = meshgrid(latV, lonV);
valG = griddata(Lo, La, Er2, latG, lonG);
contourf(latG, lonG, valG, 100, 'LineColor', 'none')
xlabel('Longitude (°)')
ylabel('Latitude (°)')
colormap(jet)
a=colorbar;
ylabel(a,'Percentage Difference %')
hold on;
bh = borders('pakistan','k-','LineWidth',2);
title('Horizontal vs Tilted Energy Difference')

