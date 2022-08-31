clear all;clc
h = 20;
w = 20;
band = 200;
show_band = 10;
r_s  = 0.1;
rho_s  = 0.1;
smooth_flag = 1;
[Omat,Lmat,Smat]=generate_M(h,w,band,r_s,rho_s,smooth_flag);
figure;
subplot(2,2,1);imshow(reshape(Lmat(:,show_band),[h,w]),[]);title('low-rank')
subplot(2,2,2);imshow(reshape(Smat(:,show_band),[h,w]),[]);title('sparse');
subplot(2,2,3);imshow(reshape(Omat(:,show_band),[h,w]),[]);title('observation');
subplot(2,2,4);imshow(reshape(Omat(:,show_band),[h,w]),[]);title('observation');

