function Image3D = load_3D_GPU(o,r,y,x,c)
%given round,y,x index of tile and channel c, this returns the full 3D image.

FileTif = o.TileFiles{r,y,x,c};
Image3D = zeros(o.TileSz,o.TileSz,o.nZ,'gpuArray');
FileID = Tiff(FileTif,'r');
rps = getTag(FileID,'RowsPerStrip');
 
for z=1:o.nZ
   setDirectory(FileID,z);
   % Go through each strip of data.
   rps = min(rps,o.TileSz);
   for j = 1:rps:o.TileSz
      row_inds = j:min(o.TileSz,j+rps-1);
      stripNum = computeStrip(FileID,j);
      Image3D(row_inds,:,z) = readEncodedStrip(FileID,stripNum);
   end
end
close(FileID);
end