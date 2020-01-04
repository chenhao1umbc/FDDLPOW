function [Dict]=FDDLOW_table2(X,trlabels,opt)

% This fucntion is designed to train the dictionary in table 2, FDDLO
% min||(D,Z,W) J(D,Z,W),  s.t.W_q^T W_q=I,q=1,2,..Q
% with J(D,Z,W)=||X-DZ||_F^2+||_1 ||Z||_1+muf(W,Z)+nug(W,Z)
% The input is  X, the training data, a matrix M by N, N data samples
%             trlabels is the labels of training data, like[1,1,1,1,2,2,3,3,3]             
%             opt are the training options with
%                 opt.K -the number of atoms in Dictionary
%                 opt.Q -the projected dimensions
%                 opt.lambda1 -lambda1 control the sparsity level
%                 opt.mu -mu the coeffecient for fisher term
%                 opt.max_iter - the maximum iteration
%                 opt.losscalc -if true then calculate loss fucntion
% The output is Dict, a struct with D,W,Z,X, Loss(the loss function value)

% initialize Dictionary
C = max(trlabels);
[M_d, N]=size(X);
Nc = N / C;
H1 = kron(eye(C),ones(Nc)/Nc);
H2 = ones(N)/N;
opt_init = opt;
opt_init.C = C; opt_init.N = N; opt_init.Nc = Nc; opt_init.M_d = M_d;
[D, Z, W, U, Loss, opt] = initdict_t2(X,trlabels,opt_init); % max_iter will change for existing dictionary


% main loop
for ii = 1:opt.max_iter   
    ii
    % update D, with U W and Z fixed
    optD = opt;
    optD.max_iter = 500;
    optD.threshold = 1e-4;
    optD.showconverge = false;
    D = DDLMD_updateD(X,optD,D,Z);
    if opt.losscalc
        Loss(1,ii) = DDLMD_Loss_mix_t2(X,trlabels,opt,W,D,Z);
    end        
    
    % update Z, with D Uand W fixed
    optZ = opt;
    optZ.max_iter = 100; % for fista
    optZ.threshold = 1e-4;
    optZ.showprogress = false; % show inside of fista
    optZ.showconverge = false; % show updateZ
    optZ.showcost= true*optZ.showprogress;
    optZ.max_Ziter = 1; % for Z update
    optZ.Zthreshold = 1e-4; 
    while 1
        Z = mix_updateZ_t2(X,trlabels,optZ, W, D, Z, U);           
        a = sum(abs(Z), 2);
        nn = sum(a ==0);
        if nn >0 
            disp('In the while loop...')
            D(:, a==0) = X(:,randi([1, N],[nn,1])); 
            Z(a==0,:) = rand(nn, N);
        else
            break;
        end
    end
    M = getM_t2(opt.K, opt.C, Nc, Z); % get M, H1, and H2 for updating W and U.
    
    if opt.losscalc
        Loss(2,ii) = DDLMD_Loss_mix_t2(X,trlabels,opt,W,D,Z);
    end
    
    % update U, with D and Z fixed.
    U = mix_updateU_t2(W, M);
    
    % update W, with D U and Z fixed
    W = mix_updateW_t2(opt, H1, H2, M, U, Z);
    if opt.losscalc
        Loss(3,ii) = DDLMD_Loss_mix_t2(X,trlabels,opt,W,D,Z);
    end    
    
    % show loss function value
    if opt.losscalc
        Loss(4,ii) = DDLMD_Loss_mix_t2(X,trlabels,opt,W,D,Z);
        Dict.Loss = Loss;
        if ii > 1            
        if abs(Loss(4, ii-1) - Loss(4, ii))/Loss(ii) < 1e-3
            break;
        end
        end
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
Dict.iter = ii;

end % end of the file