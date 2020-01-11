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
H3 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H0
H3H3t = H3*H3';  % H3 is H0 in the paper
max_eig_S = 2.1; % max(eig(M1M1t-M2M2t+1.1*eye(N)));
max_eig_H3 = 4.1667e-04; %max(eig(H3H3t));
S = 2.1*eye(N) - 2*H1 + H2;
opt_init = opt;
opt_init.C = C; opt_init.N = N; opt_init.Nc = Nc; opt_init.M_d = M_d;
[D, Z, W, U, Loss, opt] = initdict_t2(X,trlabels,opt_init); % max_iter will change for existing dictionary

optD = opt;
optD.max_iter = 500;
optD.threshold = 1e-4;
optD.showconverge = false;

optZ = opt;
optZ.max_iter = 100; % for fista
optZ.threshold = 1e-4;
optZ.showprogress = false; % show inside of fista
optZ.showconverge = false; % show updateZ
optZ.showcost= true*optZ.showprogress;

% M = getM_t2(opt.K, opt.C, Nc, Z); %debug
% loss_detail = zeros(4, 10); %debug
% r= Loss;f = r;g = f;s = r;

% main loop
for ii = 1:opt.max_iter 
%     ii
    % update D, with U W and Z fixed
    D = DDLMD_updateD(X,optD,D,Z);   
%     loss_detail(1, ii) = DDLMD_Loss_mix_t2(X,trlabels,opt,W,D,Z, M, U);
    
    % update Z, with D Uand W fixed
    while 1
        Z = mix_updateZ_t2(X, trlabels, optZ, W, D, Z, U,H3, S, H3H3t,max_eig_S, max_eig_H3);           
        a = sum(abs(Z), 2);
        nn = sum(a ==0);
        if nn >0 
            disp('In the while loop...')
            ii
            D(:, a==0) = X(:,randi([1, N],[nn,1])); 
            Z(a==0,:) = rand(nn, N);
        else
            break;
        end
    end
    M = getM_t2(opt.K, opt.C, Nc, Z); % get M, H1, and H2 for updating W and U.
%     loss_detail(2, ii) = DDLMD_Loss_mix_t2(X,trlabels,opt,W,D,Z, M, U);
    
    % update U, with D and Z fixed.
    U = mix_updateU_t2(W, M);
%     loss_detail(3, ii) = DDLMD_Loss_mix_t2(X,trlabels,opt,W,D,Z, M, U);
    
    % update W, with D U and Z fixed
    W = mix_updateW_t2(opt, S, M, U, Z); 
%     loss_detail(4, ii) = DDLMD_Loss_mix_t2(X,trlabels,opt,W,D,Z, M, U);
    
    % show loss function value
    if opt.losscalc
        wtz = W'*Z;
        fWZ= opt.mu*trace(wtz*S*wtz');
        gWZ = opt.nu* norm(W'*M - U, 'fro')^2;
        sZ = opt.lambda1*sum(abs(Z(:)));
        fid = norm(X-D*Z,'fro')^2;
        l= fid + sZ + fWZ+ gWZ;
%         g(ii) = gWZ; f(ii) = fWZ; s(ii) = sZ; r(ii) = fid;
        Loss(ii) = l ;     % DDLMD_Loss_mix_t2(X,trlabels,opt,W,D,Z, M, U)
        Dict.Loss = Loss;
        if ii > 80            
        if abs(Loss( ii-1) - Loss( ii))/Loss(ii) < 1e-4
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