%% InsertTexGraphics.m
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
 function InsertTexGraphics( searchPattern, caption, fileID )

   graphicsFiles = dir( searchPattern );

   fprintf( fileID, ['\\subsection{',caption,'}\n'] );
   fprintf( fileID, ['\\begin{figure}\n'] );
   for i = 1:length(graphicsFiles)
     fprintf( fileID, ['\\includegraphics[width=0.3\\textwidth]{',graphicsFiles(i).name,'}\n'] );
   end 
   fprintf( fileID, ['\\caption{',caption,'}\n\n'] );
   fprintf( fileID, ['\\end{figure}\n'] );




