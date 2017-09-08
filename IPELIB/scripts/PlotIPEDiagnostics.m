%% PlotIPEDiagnostics.m
%
% Author : Joseph Schoonover
%          Cooperative Institute for Research in Environmental Sciences (CIRES) 
%          National Oceanographic and Atmospheric Administration (NOAA)
%          Space Weather Prediction Center (SWPC)
%
% Copyright (2017), All Rights Reserved
%
%
% ///////////////////////////////////////////////////////////////////////////////////// %
%
%%

   % Set the min and max values for the contours, and the color bar axis for each plot

%% Electron Density at 300 km
   v = linspace( min(min(IPEStruct.ElectronDensity(:,:,43))), max(max(IPEStruct.ElectronDensity(:,:,43))), nColors );
   figure
     hold on
     contourf( IPEStruct.longitude, IPEStruct.latitude, IPEStruct.ElectronDensity(:,:,43), ...
               v, 'EdgeColor', 'None' )
     v = linspace( eDensityMin, eDensityMax, nContours ); 
     contour(  IPEStruct.longitude, IPEStruct.latitude, IPEStruct.ElectronDensity(:,:,43), ...
               v, 'EdgeColor', 'k', 'LineWidth', 2 )
  
     colormap( eDensityColorMap )
     colorbar
     caxis( [eDensityMin, eDensityMax] )
  
     xlabel( 'Geographic Longitude ( ^o E)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight','Bold' ); 
     ylabel( 'Geographic Latitude ( ^o N)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight','Bold' ); 
     title( {'Electron Density at 300 km',IPEStruct.TimeStamp}, 'FontSize', 20, 'FontName', 'Arial', 'FontWeight', 'Bold' );
     set( gca, 'LineWidth', 4, 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold' );
  
     print( ['ElectronDensity300km_',IPEStruct.TimeStamp], '-depsc2' );
   close all

   
%% TEC
%   v = linspace( min(min(IPEStruct.TEC)), max(max(IPEStruct.TEC)), nColors );
%   figure
%     hold on
%     contourf( IPEStruct.longitude, IPEStruct.latitude, IPEStruct.TEC, ...
%               v, 'EdgeColor', 'None' )
%     v = linspace( TECMin, TECMax, nContours ); 
%     contour(  IPEStruct.longitude, IPEStruct.latitude, IPEStruct.TEC, ...
%               v, 'EdgeColor', 'k', 'LineWidth', 2 )
%  
%     colormap( TECColorMap )
%     colorbar
%     caxis( [TECMin, TECMax] )
%  
%     xlabel( 'Geographic Longitude ( ^o E)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight','Bold' ); 
%     ylabel( 'Geographic Latitude ( ^o N)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight','Bold' ); 
%     title( {'TEC',IPEStruct.TimeStamp}, 'FontSize', 20, 'FontName', 'Arial', 'FontWeight', 'Bold' );
%     set( gca, 'LineWidth', 4, 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold' );
%  
%     print( ['TEC_',IPEStruct.TimeStamp], '-depsc2' );
%   close all
%
%% nmf2
   v = linspace( min(min(IPEStruct.nmf2)), max(max(IPEStruct.nmf2)), nColors );
   figure
     hold on
     contourf( IPEStruct.longitude, IPEStruct.latitude, IPEStruct.nmf2, ...
               v, 'EdgeColor', 'None' )
     v = linspace( nmf2Min, nmf2Max, nContours ); 
     contour(  IPEStruct.longitude, IPEStruct.latitude, IPEStruct.nmf2, ...
               v, 'EdgeColor', 'k', 'LineWidth', 2 )
  
     colormap( nmf2ColorMap )
     colorbar
     caxis( [nmf2Min, nmf2Max] )
  
     xlabel( 'Geographic Longitude ( ^o E)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight','Bold' ); 
     ylabel( 'Geographic Latitude ( ^o N)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight','Bold' ); 
     title( {'NmF2',IPEStruct.TimeStamp}, 'FontSize', 20, 'FontName', 'Arial', 'FontWeight', 'Bold' );
     set( gca, 'LineWidth', 4, 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold' );
  
     print( ['NmF2_',IPEStruct.TimeStamp], '-depsc2' );
   close all


%% temperature
   v = linspace( min(min(IPEStruct.temperature(:,:,43))), max(max(IPEStruct.temperature(:,:,43))), nColors );
   figure
     hold on
     contourf( IPEStruct.longitude, IPEStruct.latitude, IPEStruct.temperature(:,:,43), ...
               v, 'EdgeColor', 'None' )
     v = linspace( tempMin, tempMax, nContours ); 
     contour(  IPEStruct.longitude, IPEStruct.latitude, IPEStruct.temperature(:,:,43), ...
               v, 'EdgeColor', 'k', 'LineWidth', 2 )
  
     colormap( tempColorMap )
     colorbar
     caxis( [tempMin, tempMax] )
  
     xlabel( 'Geographic Longitude ( ^o E)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight','Bold' ); 
     ylabel( 'Geographic Latitude ( ^o N)', 'FontSize', 18, 'FontName', 'Arial', 'FontWeight','Bold' ); 
     title( {'ThermosphereTemperature',IPEStruct.TimeStamp}, 'FontSize', 20, 'FontName', 'Arial', 'FontWeight', 'Bold' );
     set( gca, 'LineWidth', 4, 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold' );
  
     print( ['ThermosphereTemperature_',IPEStruct.TimeStamp], '-depsc2' );
   close all
