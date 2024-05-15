clear all 
close all
%制御対象の伝達関数
G=tf([1],[1 1 0])
%任意の閉ループ系の極設定
Poles=[-1+i -1-i -3]
%コントローラの伝達関数導出
K = pole_assignment(G,Poles)

%閉ループ系の伝達関数
H=feedback(G,K)