%% LoadIPEStruct.m
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
 function [IPEStruct] = LoadIPEStruct( IPEStateFile )
 
  
  % Grid
  IPEStruct.z         = ncread( IPEStateFile, 'z' );
  IPEStruct.latitude  = ncread( IPEStateFile, 'latitude' );
  IPEStruct.longitude = ncread( IPEStateFile, 'longitude' );

  [IPEStruct.latitude,IPEStruct.longitude] = meshgrid( IPEStruct.latitude, IPEStruct.longitude );

  % Ions
  IPEStruct.Oplus  = ncread( IPEStateFile, 'O+' );
  IPEStruct.Hplus  = ncread( IPEStateFile, 'H+' );
  IPEStruct.Heplus = ncread( IPEStateFile, 'He+' );
  IPEStruct.Nplus  = ncread( IPEStateFile, 'N+' );
  IPEStruct.NOplus = ncread( IPEStateFile, 'NO+' );
  IPEStruct.O2plus = ncread( IPEStateFile, 'O2+' );
  IPEStruct.N2plus = ncread( IPEStateFile, 'N2+' );
  
  % Derived Quantities from Ions
  IPEStruct.ElectronDensity = ncread( IPEStateFile, 'e' );
  IPEStruct.TEC             = ncread( IPEStateFile, 'TEC' );
  IPEStruct.nmf2            = ncread( IPEStateFile, 'nmf2' );
  IPEStruct.hmf2            = ncread( IPEStateFile, 'hmf2' );


  % Neutrals
  IPEStruct.temperature = ncread( IPEStateFile, 'tn' );
  IPEStruct.u           = ncread( IPEStateFile, 'vn_zonal' );
  IPEStruct.v           = ncread( IPEStateFile, 'vn_meridional' );
  IPEStruct.w           = ncread( IPEStateFile, 'vn_vertical' );
  IPEStruct.O           = ncread( IPEStateFile, 'O' );
  IPEStruct.N2          = ncread( IPEStateFile, 'N2' );
  IPEStruct.O2          = ncread( IPEStateFile, 'O2' );


  IPEStruct.TimeStamp = IPEStateFile(end-15:end-3);
  
 

     


