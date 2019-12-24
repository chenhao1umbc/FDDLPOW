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
M2 = H1 - H2;
M1M1t = M1*M1';
M2M2t = M2*M2';
H3H3t = H3*H3';
opt.N = N;
opt.Nc = Nc;
H_bar_i = M1(1:Nc, 1:Nc);

% initialize Dictionary
[D, Z, W, U, V, Delta, Loss, opt] = initdict_t3(X, H_bar_i, H3, opt); 
% max_iter will change for existing dictionary

% main loop
for ii = 1:opt.max_iter   
    t1 = toc;
    
    % update D, with U W and Z fixed
    optD = opt;
    optD.max_iter = 500;
    optD.threshold = 1e-4;
    optD.showconverge = false;
    D = DDLMD_updateD(X,optD,D,Z);
    if opt.losscalc
        Loss(1,ii) = Loss_mix(X, H_bar_i, H3, M1, M2, opt,W,D,Z,U,V,Delta);
    end
        
    % update Z, with D Uand W fixed
    optZ = opt;
    optZ.M1M1t = M1M1t;
    optZ.M2M2t = M2M2t;
    optZ.H3H3t = H3H3t;
    optZ.M1 = M1;
    optZ.M2 = M2;
    
    optZ.max_iter = 500; % for fista
    optZ.threshold = 1e-6;
    optZ.showprogress = false; % show inside of fista
    optZ.showconverge = false; % show updateZ
    optZ.showcost= true*optZ.showprogress;
    optZ.max_Ziter = 20; % for Z update
    optZ.Zthreshold = 1e-6;   
    while 1
        Z = mix_updateZ(X,H_bar_i, H3, optZ, W, D, Z, U, V, Delta); 
        a = sum(abs(Z), 2);
        nn = sum(a ==0);
        if  nn >0 
%             disp('In the while loop...')
            D(:, a==0) = X(:,randi([1, N],[nn,1])); 
            Z = randn(size(Z));
        else
            break;
        end
    end
    
    if 0.3 == ii/opt.max_iter
        sparsity = mean(sum(Z ~= 0))/opt.K       % avg number of nonzero elements in cols of Z
        if sparsity > 0.9 || sparsity < 0.05
            fprintf('30 percent iters too sparse or non-sparse\n')
            break;            
        end
    end
    if opt.losscalc
        Loss(2,ii) = Loss_mix(X, H_bar_i, H3, M1, M2, opt,W,D,Z,U,V,Delta);
    end
    
    % update W
    W = mix_updateW(opt,H_bar_i, M1, M2, H3, Delta, U, V,  Z);    
    if opt.losscalc
        Loss(3,ii) = Loss_mix(X, H_bar_i, H3, M1, M2, opt,W,D,Z,U,V,Delta);
    end     
    
    % update U, with D and Z fixed.    
    U = mix_updateU(W, Z, H3);

    % updtae V
    V = mix_updateV(H_bar_i, Z, W, Delta, opt);

    % update Delta   
    Delta = mix_updateDelta(H_bar_i, Z, W, V, opt);

% a = cell(1, C);
% b = cell(1, C);
% dist = 0;
% for i = 1:C    
%     a{i} = H_bar_i*(W'*Z(:, 1+ Nc*(i-1): Nc*i))';
%     b{i} = Delta(i)*V{i};
%     dist = dist + norm(a{i} - b{i}, 'fro')^2; % perclass whitening term
% end
% dist
% normW = norm(W, 'fro')

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