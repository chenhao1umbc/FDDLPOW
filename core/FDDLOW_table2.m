function [Dict]=FDDLOW_table2(X,trlabels,opt)
% This fucntion is designed to train the dictionary in table 2, FDDLO
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
[D, Z, W, U, V, Delta, Loss, opt] = initdict(X,trlabels,opt); % max_iter will change for existing dictionary

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
        Loss(1,ii) = DDLMD_Loss_mix(X,trlabels,opt,W,D,Z,U,V,Delta);
    end
        
    % update Z, with D Uand W fixed
    optZ = opt;
    optZ.max_iter = 300; % for fista
    optZ.threshold = 1e-6;
    optZ.showprogress = false; % show inside of fista
    optZ.showconverge = false; % show updateZ
    optZ.showcost= true*optZ.showprogress;
    optZ.max_Ziter = 10; % for Z update
    optZ.Zthreshold = 1e-6;        
    Z = mix_updateZ(X,trlabels,optZ, W, D, Z, U, V, Delta); 
    [H1, H2, H3] = getMH1H2H3(trlabels, Z); % get M, H1, and H2 for updating W and U.
    sparsity = mean(sum(Z ~= 0))       % avg number of nonzero elements in cols of Z
    if opt.losscalc
        Loss(2,ii) = DDLMD_Loss_mix(X,trlabels,opt,W,D,Z,U,V,Delta);
    end
    
    % update U, with D and Z fixed.
    U = mix_updateU(W, Z, Delta, H3);

    % update Delta
    Delta = sum(sum(U.*(W'*Z*H3)))/norm(U,'fro')^2;    
    
    % updtae V
    V = mix_updateV(W, Z, H1);
    
    % update W, with D Uand Z fixed
    W = mix_updateW(opt, H1, H2, H3, Delta, U, V,  Z);    
    if opt.losscalc
        Loss(3,ii) = DDLMD_Loss_mix(X,trlabels,opt,W,D,Z,U,V,Delta);
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
    if opt.savedict
        if mod(ii,60) == 0
            Dict.D = D;
            Dict.W = W;
            Dict.Z = Z;
            Dict.U = U;
            Dict.V = V;
            Dict.Delta = Delta;
            Dict_mix = Dict;
            opts = opt;
            save([opts.mixnm(1:end-4),'_',num2str(ii)],'Dict_mix','opts')
        end
    end
    
fprintf('one iter time: %6.4f \n',toc-t1)
end

Dict.D = D;
Dict.W = W;
Dict.Z = Z;
Dict.U = U;
Dict.V = V;
Dict.Delta = Delta;

end % end of the file