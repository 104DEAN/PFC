clear all 
close all

G=tf([1],[1 1 0])
Poles=[-1+i -1-i -3]
K = pole_assignment(G,Poles)