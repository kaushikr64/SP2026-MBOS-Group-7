function [F, DF] = CR3BP_PerOrb_SF(X_guess, mu, section)

% Unpack terms for poincare section
X_sec = section.X_sec;
n_sec = section.n_sec/norm(section.n_sec);
null_sec = null(n_sec');

% Integrate state to Poincare section return
[X_return, DP_return, STM_return, t_return] = CR3BP_Prop2Poincare(X_guess, mu, section, [0,20]);


F = [null_sec'*(X_return-X_guess);... % Residual between return state and guess on the section)
    n_sec'*(X_guess-X_sec)]; % Residual of the guess from the section

DF = [null_sec'*(DP_return - eye(6)); n_sec'];
end