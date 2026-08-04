*
* Grid for a rounded expanding corner
*
* $gridName = "expansionCorner1.hdf"; $factor=1;
$gridName = "expansionCorner2.hdf"; $factor=2;
* $gridName = "expansionCorner3.hdf"; $factor=4;
* $gridName = "expansionCorner4.hdf"; $factor=8;
*
* 
create mappings
  smoothedPolygon
    vertices
    3
    -.5   0.     
     0.   0.     
     0.4  -.2  
*     3
*     -1. 0. 
*     0. 0.
*     .7 -.7 
    n-dist
    fixed normal distance
      .3 .5 1. 
    sharpness
      10
      10
      10
    t-stretch
      0 50
      .1 15.
      0 50
    n-stretch
      1. 3. 0
    correct corners
    lines
      *        75 25 
      $nx = int(40*$factor+1);  $ny=int(15*$factor+1);
      $nx $ny
    mappingName
    corner
    exit
exit
*
generate an overlapping grid
  corner
  done
  change parameters
    ghost points
      all
      2 2 2 2 2 2
  exit
  compute overlap
exit
*
save an overlapping grid
  $gridName
  expansionCorner
exit





















*
* Grid for a rounded expanding corner
*
create mappings
  smoothedPolygon
    vertices
    3
    -.5   0.     
     0.   0.     
     0.4  -.2  
*     3
*     -1. 0. 
*     0. 0.
*     .7 -.7 
    n-dist
    fixed normal distance
      .3 .5 1. 
    sharpness
    10
    10
    10
    t-stretch
    0 50
    .1 20.
    0 50
    n-stretch
    1. 2. 0
    correct corners
    lines
       37 13  * 75 25 41
    mappingName
    corner
    exit
exit
*
generate an overlapping grid
  corner
  done
  change parameters
    ghost points
      all
      2 2 2 2 2 2
  exit
  compute overlap
exit
*
save an overlapping grid
  expansionCorner.hdf
  expansionCorner
exit
