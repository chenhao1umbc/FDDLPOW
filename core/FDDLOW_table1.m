function [Dict]=FDDLOW_table1(X,trlabels,opt)
% This fucntion is designed to train the dictionary in (8)
% min||(D,Z,W) J(D,Z,W),  s.t.W_q^T W_q=I,q=1,2,..Q
% with J(D,Z,W)=||X-DZ||_F^2+||_1 ||Z||_1+muf(W,Z)
% f(W,Z)=tr(W^T{S_W}W) + tr(W^T{S_B}W) + ||{W^T}Z||_F^2
% The input is  X, the training data, a matrix M by N, N data samples
%             trlabels is the labels of training data, like[1,1,1,1,2,2,3,3,3]             
%             opt are the training options with
%                 opt.K -the number of atoms in Dictionary
%                 opt.Q -the projected dimensions
%                 opt.lambda1 -lambda1 control the sparsity level
%                 opt.mu -mu the coeffecient for fisher term
%                 opt.max_iter - the maximum iteration
%                 opt.losscalc -if true then calculate loss fucntion
% The output is Dict, a struct with D,W,Z, Loss(the loss function value)

% initialize Dictionary
[D, Z, W, Loss, opt]=initdict(X,trlabels,opt); % max_iter will change for existing dictionary

% main loop
for ii=1:opt.max_iter  
    ii
    %     tic
        % update D, with W and Z fixed
        optD=opt;
        optD.max_iter=500;
        optD.threshold=5e-5;
        optD.showconverge=false;
        D=DDLMD_updateD(X,optD,D,Z);

    while 1
        % update Z, with D and W fixed
        optZ=opt;
        optZ.max_iter=500;  % for fista iters
        optZ.threshold=1e-6;
        optZ.showprogress = false; % show inside of fista
        optZ.showconverge = false; % show updateZ
        optZ.showcost= true*optZ.showprogress;         
        Z=DDLMD_updateZ(X,trlabels,optZ,W,D,Z);
        
        a = sum(abs(Z), 2);
        if sum(a ==0) >0 
            disp('In the while loop...')
            D(:, a==0) = X(:,randi([1, 800],[sum(a ==0),1])); 
            Z = randn(size(Z));
        else
            break;
        end
    end
    if 0.3 == ii/opt.max_iter
        sparsity=mean(sum(Z ~= 0))/opt.K;   % avg number of nonzero elements in cols of Z
        if sparsity > 0.95 || sparsity < 0.05
            fprintf('30 percent iters too sparse or non-sparse\n')
            break;            
        end
    end
    % update W, with D and Z fixed
    optW=opt;
    optW.ploteig=false;
    W=DDLMD_updateW(trlabels,optW,Z);    
      

    % show loss function value
    if opt.losscalc
        Loss(ii)=DDLMD_Loss(X,trlabels,opt,W,D,Z);
        Dict.Loss=Loss;
        if ii > 1            
        if abs(Loss(ii-1) - Loss(ii)) < 1e-5
            break;
        end
        end
    end
    if opt.showconverge
        figure(400);
        plot(Loss(:));
        title('objective function value');
        xlabel('Iterations');
        pause(.1);
    end
    
    if opt.savedict*0
        if mod(ii,40)==0
            Dict.D=D;
            Dict.W=W;
            Dict.Z=Z;
            opts=opt;
            save([opts.Dictnm(1:end-4),'_',num2str(ii)],'Dict','opts')
        end
    end
%     toc
end

Dict.D=D;
Dict.W=W;
Dict.Z=Z;
Dict.iter = ii;

end % end of the function file