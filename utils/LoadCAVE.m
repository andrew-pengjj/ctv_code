function [Z]  =  LoadCAVE(image,Dir)

% Oct 2021
% written by Jiangjun Peng

fpath         =   fullfile( fullfile(Dir, image, image), '*.png');
im_dir        =   dir(fpath);
bands         =   length(im_dir);
filestem      =   char(strcat(Dir, '/',image, '/', image, '/', image));
for band = 1:bands
    filesffx    =  '.png';    
    if band < 10
        prefix  =  '_0';
    else
        prefix  =  '_';
    end    
    number      =   strcat( prefix, int2str(band) );
    filename    =   strcat( filestem, number, filesffx );
    Z           =   double(imread(filename));
    if length(size(Z))==3
        Z = mean(Z,3);
    end
    if band ==1
        sz = size(Z);
    end
    maxv          =   max(Z(:));
    minv          =   min(Z(:));
    s_Z(band, :)=   (Z(:)-minv)/(maxv-minv); 
end
s_Z      =  s_Z';
Z        = reshape(s_Z,[sz(1),sz(2),bands]);
end