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
   fprintf( fileID, ['\\begin{center}\n'] );
   for i = 1:length(graphicsFiles)
     fprintf( fileID, ['\\includegraphics[width=0.3\\textwidth]{',graphicsFiles(i).name,'}\n'] );
     if( mod( i, 3 ) == 0 )
       fprintf( fileID, ['\\end{center}\n'] );
       fprintf( fileID, ['\\end{figure}\n\n'] );
       fprintf( fileID, ['\\begin{figure}\n'] );
       fprintf( fileID, ['\\begin{center}\n'] );
     end
   end 

   if( mod( i, 3 ) ~= 0 )
     fprintf( fileID, ['\\end{figure}\n'] );
   end




