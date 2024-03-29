function [Dict]=FDDLOW_table3(X, opt)

% This fucntion is designed to train the dictionary in table 3, FDDLOW
%
%
% The input is  X, the training data, a matrix M by N, N data samples
%             trlabels is the labels of training data, like[1,1,1,1,2,2,3,3,3]             
%             opt are the training options with
%                 opt.K -the number of atoms in Dictionary
%                 opt.Q -the projected dimensions
%                 opt.lambda1 -lambda1 control the sparsity level
%                 opt.mu -mu the coeffecient for fisher term
%                 opt.max_iter - the maximum iteration
%                 opt.losscalc -if true then calculate loss fucntion
% The output is Dict, a struct with D,W,Z,X,Delta,U,V, Loss(the loss function value)

C = opt.C;
[~, N]=size(X);
Nc =N/C;
[H1, H2, H3] = getMH1H2H3(N, Nc, C);
M1 = eye(N) - H1;
max_eig_S = 2.1;
max_eig_H3 = 4.1667e-04;
sum_max_eig_Hhat_c = 6;
H3H3t = H3*H3';
opt.N = N;
opt.Nc = Nc;
H_bar_i = M1(1:Nc, 1:Nc);  
S = 2.1*eye(N) - 2*H1 + H2;
% H_hat = cell(opt.C);
% temp = zeros(N, Nc);
% for ii = 1:opt.C
%     temp(1+ Nc*(ii-1): Nc*ii,:) = H_bar_i;
%     H_hat{ii} = temp;
%     temp = temp -temp;
% end  

% initialize Dictionary
[D, Z, W, U, V, Delta, Loss, opt] = initdict_t3(X, H_bar_i, H3, opt); 
% max_iter will change for existing dictionary

optD = opt;
optD.max_iter = 500;
optD.threshold = 1e-4;
optD.showconverge = false;

optZ = opt;
optZ.S = S;
optZ.max_eig_S = max_eig_S;
optZ.max_eig_H3 = max_eig_H3;
optZ.H_bar_iSq = H_bar_i^2;
optZ.H3H3t = H3H3t;
optZ.sum_max_eig_Hhat_c = sum_max_eig_Hhat_c;
optZ.max_iter = 200; % for fista
optZ.threshold = 1e-4;
optZ.showprogress = false; % show inside of fista
optZ.showconverge = false; % show updateZ
optZ.showcost= true*optZ.showprogress;

% main loop
for ii = 1:opt.max_iter       
    % update D, with U W and Z fixed
    D = DDLMD_updateD(X,optD,D,Z);
       
    % update Z, with D Uand W fixed 
    while 1
        Z = mix_updateZ(X,H_bar_i, H3, optZ, W, D, Z, U, V, Delta); 
        a = sum(abs(Z), 2);
        nn = sum(a ==0);
        if  nn >0 
            ii
            disp('In the while loop...')
            D(:, a==0) = X(:,randi([1, N],[nn,1])); 
            Z(a==0,:) = rand(nn, N);
        else
            break;
        end
    end
    
    % update W
    W = mix_updateW(opt,H_bar_i,S, H3, Delta, U, V,  Z);      

    % update U, with D and Z fixed.    
    U = mix_updateU(W, Z, H3);

    % updtae V
    [V, WtZHhatc] = mix_updateV(H_bar_i, Z, W, Delta, opt);

    % update Delta   
    Delta = mix_updateDelta(WtZHhatc, V, opt);
    
    Loss(ii) = Loss_mix(X, H_bar_i, H3,S, opt,W,D,Z,U,V,Delta);
    Dict.Loss = Loss;
    if ii > 120            
    if abs(Loss( ii-1) - Loss( ii))/Loss(ii) < 1e-4
        break;
    end
    end
    % show loss function value
    if opt.losscalc
        Dict.Loss = Loss;
        if opt.showconverge 
            figure(800);            
            semilogy(vec(Loss(:,1:ii)));        
            title('Cost function for mixture case');
            xlabel({'Iterations';'--from DDLMD\_mix.m'});
            pause(.1);
        end
    end
    
% fprintf('one iter time: %6.4f \n',toc-t1)
end

Dict.D = D;
Dict.W = W;
Dict.Z = Z;
Dict.U = U;
Dict.V = V;
Dict.Delta = Delta;
Dict.iter = ii;

end % end of the file