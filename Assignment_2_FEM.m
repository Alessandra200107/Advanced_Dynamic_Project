%% ASSIGNMENT 2 - FEM

clear all
close all
clc

% load input file and assemble structure
[file_i,xy,nnod,sizee,idb,ndof,incidenze,l,gamma,m,EA,EJ,posiz,nbeam,pr]=loadstructure;

% draw structure
dis_stru(posiz,l,gamma,xy,pr,idb,ndof);
