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
N = size(X, 2);
% main loop
for ii=1:opt.max_iter  
    ii
    % update D, with W and Z fixed
%     disp('Updating D ...');
    optD=opt;
    optD.max_iter=500;
    optD.threshold=5e-5;
    optD.showconverge=false;
    D_old = D;
    D=DDLMD_updateD(X,optD,D,Z);

    % update Z, with D and W fixed
%     disp('Updating Z ...');
    optZ=opt;
    optZ.max_iter=100;  % for fista iters
    optZ.threshold=1e-4;    
    optZ.showprogress = false; % show updateZ
    optZ.showcost= optZ.showprogress;  
    Z_old = Z;
    while true
        Z = DDLMD_updateZ(X,trlabels,optZ,W,D,Z);    
        a = sum(abs(Z),2);
        nn = sum(a ==0);
        if nn > 0 
            disp('Expunging unused atoms ...')
            find(a == 0)
            D(:,a==0) = X(:,randi([1, N],[nn,1])); 
            Z(a==0,:) = rand(nn,N);
        else
            break;
        end
    end
    
    % update W, with D and Z fixed
    optW=opt;
    optW.ploteig=false;
    W_old = W;
    W=DDLMD_updateW(trlabels,optW,Z);    
      

    % show loss function value
    if opt.losscalc
        Loss(ii)=DDLMD_Loss(X,trlabels,opt,W,D,Z);
        Dict.Loss=Loss;
        if ii > 1
            if abs(Loss(ii-1) - Loss(ii))/Loss(ii) < 1e-4
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
    
    % stopping criterion
    deltaD = max(abs(D(:) - D_old(:)));
    deltaZ = max(abs(Z(:) - Z_old(:)));
    deltaW = max(abs(W(:) - W_old(:)));
    delta = max([deltaD,deltaZ,deltaW]);
    if delta < opt.th
        break;
    end
    
%     if opt.savedict
%         if mod(ii,40)==0
%             Dict.D=D;
%             Dict.W=W;
%             Dict.Z=Z;
%             opts=opt;
%             save([opts.Dictnm(1:end-4),'_',num2str(ii)],'Dict','opts')
%         end
%     end
end

Dict.D=D;
Dict.W=W;
Dict.Z=Z;
Dict.iter = ii;

end % end of the function file