%% StartTexFile.m
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
  
 texFid = fopen( [PlotDir,'Report.tex'], 'w' );

 fprintf( texFid, '\\documentclass[12pt,a4paper]{article}\n' );
 fprintf( texFid, '\\usepackage[utf8]{inputenc}\n' );
 fprintf( texFid, '\\usepackage{amsmath}\n' );
 fprintf( texFid, '\\usepackage{amsfonts}\n' );
 fprintf( texFid, '\\usepackage{amssymb}\n' );
 fprintf( texFid, '\\usepackage{graphicx}\n' );
 fprintf( texFid, '\\usepackage[margin=3cm]{geometry}\n');
 fprintf( texFid, ['\\title{',model,' \\\\ ',tag,' Updates and assessments}\n']);
 fprintf( texFid, '\\begin{document}\n');
 fprintf( texFid, '\\maketitle\n');
 fprintf( texFid, '\\section{Summary of updates}\n\n');
 fprintf( texFid, '\\section{Model Output}\n\n');

     


