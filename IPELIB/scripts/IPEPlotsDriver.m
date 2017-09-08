%% IPEPlotsDriver.m
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
%% User Input

 function IPEPlotsDriver( model, tag, IPEDir, PlotDir )


  fprintf( ['Model : ', model, '\n'])
  fprintf( ['Tag   : ', tag, '\n'])
  fprintf( ['NetCDF: ', IPEDir, '\n'])
  fprintf( ['Plot  : ', PlotDir, '\n'])

  % Plotting Settings
  nColors = 200;
  nContours = 7;


  % Electron Density at 300 km
  eDensityMin = 3.0*10^8;
  eDensityMax = 3.0*10^12;
  eDensityColorMap = reds( nColors );

  % TEC
  TECMin = 3.0*10^15;
  TECMax = 8.0*10^17;
  TECColorMap = reds( nColors );

  % nmf2
  nmf2Min = 2.0*10^10;
  nmf2Max = 4.0*10^12;
  nmf2ColorMap = reds( nColors );

  % thermosphere temperature
  tempMin = 750;
  tempMax = 1100;
  tempColorMap = reds( nColors );



%% Setup

  IPEfiles = dir( [IPEDir, 'IPE_State.*.nc'] );

  StartTexFile

%% Make the Plots

  for i = 1:length(IPEfiles);

     IPEStruct = LoadIPEStruct( [IPEDir,IPEfiles(i).name] );
     PlotIPEDiagnostics 

  end

  % Insert Electron Density at 300km
  InsertTexGraphics( [PlotDir,'ElectronDensity300km_*.eps'],...
                     'Electron Density at 300 km', texFid )
  % Insert TEC
  InsertTexGraphics( [PlotDir,'TEC_*.eps'],...
                     'Total Electron Content', texFid )

  % Insert NmF2
  InsertTexGraphics( [PlotDir,'NmF2_*.eps'],...
                     'NmF2', texFid )

  % Insert Thermosphere temperature at 300 km
  InsertTexGraphics( [PlotDir,'ThermosphereTemperature_*.eps'],...
                     'Thermosphere Temperature at 300 km', texFid )

  FinalizeTexFile

     
  

