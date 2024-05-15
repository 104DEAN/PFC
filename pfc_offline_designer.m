function [k_0, nu_x] = pfc_offline_designer(sys, T_CLRT, n_B, n_h)
% --- PFC offline designer ---
% ver:0.0.3 (2022-04-01, MATLAB R2020b)
% author:KASAHARA Junpei
%
% PFCによる制御を行う際に必要となる値を，制御対象の次数に合わせて計算する．
%
% 返り値はセットポイントとプラント出力の差に掛かる係数k_0と，
% 内部モデルの状態変数に掛かる係数ベクトルnu_x．
% 制御計算に使用するのはnu_xの転置だが，返り値は転置していないことに注意．
%
% [k_0, nu_x] = pfc_offline_designer(sys, T_CLRT, n_B, n_h);
% --- Arguments ---
% sys   : ss               離散時間状態空間モデル
% T_CLRT: {mustBeNumeric}  閉ループ応答時間[s]
% n_B   : {mustBeNumeric}  基底関数の個数
% n_h   : {mustBeNumeric}  一致点の個数
% ---- Returns ----
% k_0   : double
% nu_x  : double (1, Nx)   ※Nxは状態変数の個数
%
% 【参考文献】
%  佐藤 俊之，「シンプルなモデル予測制御"Predictive Functional Control"」，
%  設計工学，Vol.52，No.10，2017．
    arguments
        sys    ss               % 離散時間状態空間モデル
        T_CLRT {mustBeNumeric}  % 閉ループ応答時間[s]
        n_B    {mustBeNumeric}  % 基底関数の個数
        n_h    {mustBeNumeric}  % 一致点の個数
    end
    
    % 計算条件確認
    if sys.Ts == 0 
        % 連続時間モデルの場合はエラーで終了．事前に離散化しておくこと．
        error('The input ss object is a continuous-time model (Ts=0). Please enter the discrete-time model.');
    elseif n_B == 0
        error('n_B == 0');
    elseif n_h == 0
        error('n_h == 0');
    end
    
    % 一致点のサンプル時刻をまとめた配列 (1, n_h)
    h = floor( T_CLRT ./ (sys.Ts *  (n_h:-1:1).') );

    % 参照軌道の減衰率
    % （分子係数の3は厳密にはTsとT_CLRTの関数として表されるが，T_CLRT≒0で
    % 　特異点に落ち込む以外は殆ど3に近い値を取るので，3に固定して計算する）
    alpha = exp(-3 * sys.Ts / T_CLRT);
    
    % nuとnu_xの計算に使用する行列を計算
    tmp0 = zeros(n_B, n_h);
    tmp1 = zeros(n_B, n_B);
    tmp2 = zeros(n_h, length(sys.A));
    for j_cnt = 1:n_h
        y_B = calc_y_B(sys, n_B, h(j_cnt));
        tmp0(:, j_cnt) = y_B;
        tmp1 = tmp1 + y_B * y_B.';
        tmp2(j_cnt, :) = sys.C * sys.A^h(j_cnt) - sys.C;
    end
    tmp1 = inv(tmp1);
    nu = tmp0.' * tmp1(:, 1);

    % 定数k_0を計算
    k_0 = nu.' * (1 - alpha.^h);
    
    % 定数ベクトルν_xを計算
    nu_x = -tmp2.' * nu;
end

% ===================== Local Function ====================== %

% 一致点における各基底関数に対するモデル出力をまとめた
% ベクトル（参考文献の式(13)）を計算する．
%
% --- Arguments ---
% sys : ss               離散時間状態空間モデル
% n_B : {mustBeNumeric}  基底関数の個数
% h_j : {mustBeNumeric}  一致点のサンプル時刻
% ---- Return ----
% y_B : double (n_B, 1)
%--------------------------------------------------------------------------
function y_B = calc_y_B(sys, n_B, h_j)
    arguments
       sys    ss               % 離散時間状態空間モデル
       n_B    {mustBeNumeric}  % 基底関数の個数
       h_j    {mustBeNumeric}  % 一致点の個数
    end

    t = (0:h_j) * sys.Ts;  % 時間軸ベクトル
    y_B = zeros(n_B, 1);
    for l = 1:n_B
        % --- 基底関数入力に対する応答を計算 --- %
        % 参考文献の式(21)をそのまま実装するよりも5倍ほど高速に計算できる（計算結果は同じ）．
        u = (0:h_j).^(l-1);  % 基底関数入力
        y = lsim(sys, u, t, 'zoh');
        % ------------------------------------ %
        
        y_B(l) = y(h_j + 1);  % 一致点における出力だけ取り出す
    end
end
% ============================================================ %