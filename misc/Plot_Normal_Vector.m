%This script is used to visulize normal vectors for fault nodes generated by EQdyna2d_2.0.0.
% 09082020, DL.
% Input: vert.txt, nsmpnv.txt, nsmp.txt
%  x1.txt, x2.txt, x3.txt ... (or other names for files storing fault geometry contolling points)
clear all; close all;

nft = [295,195,1769];
nmax = 1769;
path_mesh = './XX/' %please define the path to mesh files: vert.txt, etc. 
path_cp = './XX/' %please define the path to controlling points files: x1.txt, etc. 

x1 = load(strcat(path_cp,'/x1_1.txt'));
x1 = load(strcat(path_cp,'/x2_1.txt'));
x1 = load(strcat(path_cp,'/x3_1.txt'));

figure(1)

plot(x1(:,1),x1(:,2),'m');
plot(x2(:,1),x2(:,2),'c');
plot(x3(:,1),x3(:,2),'r');

nsmpnv = load(strcat(path_mesh,'nsmpnv.txt'));
vec = load(strcat(path_mesh,'vert.txt')); vec = vec/1000;
nsmp = load(strcat(path_mesh,'nsmp.txt'));

for i = 1: nft(1)
    quiver(vec(nsmp(i,1),1), vec(nsmp(i,1),2), nsmpnv(i,1), nsmpnv(i,2)); hold on;
end

for i = 1+nmax: nmax + nft(2)
    quiver(vec(nsmp(i,1),1), vec(nsmp(i,1),2), nsmpnv(i,1), nsmpnv(i,2)); hold on;
end

for i = 1+2*nmax: 2*nmax + nft(3)
    quiver(vec(nsmp(i,1),1), vec(nsmp(i,1),2), nsmpnv(i,1), nsmpnv(i,2)); hold on;
end
title('Fault geometry and normal vectors of split nodes');
xlabel('km');
ylabel('km');
axis equal;