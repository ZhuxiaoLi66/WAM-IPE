
function [map] = reds( nColors )


   map = ones(nColors,3);

   map(:,2) = linspace( 1, 0, nColors );
   map(:,3) = linspace( 1, 0, nColors );
