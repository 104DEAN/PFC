function [k_0, nu_x] = pfc_offline_designer(sys, T_CLRT, n_B, n_h)
% --- PFC offline designer ---
% ver:0.0.3 (2022-04-01, MATLAB R2020b)
% author:KASAHARA Junpei
%
% PFC�ɂ�鐧����s���ۂɕK�v�ƂȂ�l���C����Ώۂ̎����ɍ��킹�Čv�Z����D
%
% �Ԃ�l�̓Z�b�g�|�C���g�ƃv�����g�o�͂̍��Ɋ|����W��k_0�ƁC
% �������f���̏�ԕϐ��Ɋ|����W���x�N�g��nu_x�D
% ����v�Z�Ɏg�p����̂�nu_x�̓]�u�����C�Ԃ�l�͓]�u���Ă��Ȃ����Ƃɒ��ӁD
%
% [k_0, nu_x] = pfc_offline_designer(sys, T_CLRT, n_B, n_h);
% --- Arguments ---
% sys   : ss               ���U���ԏ�ԋ�ԃ��f��
% T_CLRT: {mustBeNumeric}  ���[�v��������[s]
% n_B   : {mustBeNumeric}  ���֐��̌�
% n_h   : {mustBeNumeric}  ��v�_�̌�
% ---- Returns ----
% k_0   : double
% nu_x  : double (1, Nx)   ��Nx�͏�ԕϐ��̌�
%
% �y�Q�l�����z
%  ���� �r�V�C�u�V���v���ȃ��f���\������"Predictive Functional Control"�v�C
%  �݌v�H�w�CVol.52�CNo.10�C2017�D
    arguments
        sys    ss               % ���U���ԏ�ԋ�ԃ��f��
        T_CLRT {mustBeNumeric}  % ���[�v��������[s]
        n_B    {mustBeNumeric}  % ���֐��̌�
        n_h    {mustBeNumeric}  % ��v�_�̌�
    end
    
    % �v�Z�����m�F
    if sys.Ts == 0 
        % �A�����ԃ��f���̏ꍇ�̓G���[�ŏI���D���O�ɗ��U�����Ă������ƁD
        error('The input ss object is a continuous-time model (Ts=0). Please enter the discrete-time model.');
    elseif n_B == 0
        error('n_B == 0');
    elseif n_h == 0
        error('n_h == 0');
    end
    
    % ��v�_�̃T���v���������܂Ƃ߂��z�� (1, n_h)
    h = floor( T_CLRT ./ (sys.Ts *  (n_h:-1:1).') );

    % �Q�ƋO���̌�����
    % �i���q�W����3�͌����ɂ�Ts��T_CLRT�̊֐��Ƃ��ĕ\����邪�CT_CLRT��0��
    % �@���ٓ_�ɗ������ވȊO�͖w��3�ɋ߂��l�����̂ŁC3�ɌŒ肵�Čv�Z����j
    alpha = exp(-3 * sys.Ts / T_CLRT);
    
    % nu��nu_x�̌v�Z�Ɏg�p����s����v�Z
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

    % �萔k_0���v�Z
    k_0 = nu.' * (1 - alpha.^h);
    
    % �萔�x�N�g����_x���v�Z
    nu_x = -tmp2.' * nu;
end

% ===================== Local Function ====================== %

% ��v�_�ɂ�����e���֐��ɑ΂��郂�f���o�͂��܂Ƃ߂�
% �x�N�g���i�Q�l�����̎�(13)�j���v�Z����D
%
% --- Arguments ---
% sys : ss               ���U���ԏ�ԋ�ԃ��f��
% n_B : {mustBeNumeric}  ���֐��̌�
% h_j : {mustBeNumeric}  ��v�_�̃T���v������
% ---- Return ----
% y_B : double (n_B, 1)
%--------------------------------------------------------------------------
function y_B = calc_y_B(sys, n_B, h_j)
    arguments
       sys    ss               % ���U���ԏ�ԋ�ԃ��f��
       n_B    {mustBeNumeric}  % ���֐��̌�
       h_j    {mustBeNumeric}  % ��v�_�̌�
    end

    t = (0:h_j) * sys.Ts;  % ���Ԏ��x�N�g��
    y_B = zeros(n_B, 1);
    for l = 1:n_B
        % --- ���֐����͂ɑ΂��鉞�����v�Z --- %
        % �Q�l�����̎�(21)�����̂܂܎����������5�{�قǍ����Ɍv�Z�ł���i�v�Z���ʂ͓����j�D
        u = (0:h_j).^(l-1);  % ���֐�����
        y = lsim(sys, u, t, 'zoh');
        % ------------------------------------ %
        
        y_B(l) = y(h_j + 1);  % ��v�_�ɂ�����o�͂������o��
    end
end
% ============================================================ %