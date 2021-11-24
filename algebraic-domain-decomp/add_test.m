clear;
clc;
close all;

%% Test 1 (bucky ball)
% B = bucky;
% figure()
% spy(B)
% title('SPSD Matrix Before CM Reordering')
% Bcm = B(symrcm(B),symrcm(B));
% Bcm = Bcm + 100*eye(size(Bcm));
% figure()
% spy(Bcm)
% title('SPSD Matrix After CM Reordering')
% A = abs(Bcm);


%% Test 2 (Semi Random Matrix)
% load west0479.mat
% B = west0479;
% S = B' * B + 50*speye(size(B));
% S = abs(S);
% pct = 100 / numel(S);
% 
% figure()
% spy(S)
% title('SPSD Matrix Before CM Reordering')
% A=(S(symrcm(S),symrcm(S)));
% figure()
% spy(A)
% title('SPSD Matrix After CM Reordering')

%% Test 3 (airfoil)
load airfoil;
% Scaling x and y
x = pow2(x,-32); 
y = pow2(y,-32);

% Forming the sparse adjacency matrix and making sure it is positive definite
n = max(max(i),max(j));
A = sparse(i,j,-1,n,n);
A = A + A';
d = abs(sum(A)) + 1;
A = A + diag(sparse(d));
spy(A)
title('SPSD Matrix Before CM Reordering')
B=A;
Bcm=(A(symrcm(A),symrcm(A)));
figure()
spy(Bcm)
title('SPSD Matrix After CM Reordering')

% Number of tests to run
q = 20;

% Allocate timing slots
fulltimes = zeros([1,q]);
blocktimes = zeros([1,q]);

for z=1:q
    % Divide Into Symmetric Diagonal Blocks
    A = abs(A);
    n = length(A);
    nblocks = 48;
    nvert = nblocks+1;
    vert = zeros([1,nvert]);
    vert(1) = 1;
    vert(nvert)=n;
    for k=1:nvert-2
        vert(k+1) = round(((k)*n/nblocks));
    end
    vert;
    % Find first block A1
    bands = find(A(vert(2),:));
    rbandwidth = bands(length(bands)) - vert(2);
    eval('A1 = A(1:vert(2)+rbandwidth,1:vert(2)+rbandwidth);');

    % Find kth block Ak
    for k=2:nblocks-1
        bands = find(A(vert(k+1),:));
        rbandwidth = bands(length(bands)) - vert(k+1);
        eval(strcat('A', num2str(k), '= A(vert(k):vert(k+1)+rbandwidth,vert(k):vert(k+1)+rbandwidth);'));
    end

    % Find nth block An
    tempvertex=vert(length(vert)-1);
    eval(strcat('A', num2str(k+1), '= A(tempvertex:n,tempvertex:n);'));
    
    title('SPSD Matrix after CM Reordering')

    % Generate a random f vector
    f = rand(length(A),1);
    
    % Time PCG on the full system
    tic;
    eval(strcat('L = ichol(A,struct(''type'',''ict'',''droptol'',1e-5));'));
    eval('pcg(A,f,10e-5,1000000,L,L'');');
    fulltimeelapsed=toc;
    blocktimeelapsed=zeros([1,nblocks]);
    for k=1:nblocks
        
        % Time the kth block
        tic;
            % FULL CHOLESKY
                eval(strcat('L = chol(A',num2str(k),');'));
            % INCOMPLETE CHOLESKY WITH DROPPING
                % eval(strcat('L = ichol(A',num2str(k),',struct(''type'',''ict'',''droptol'',1e-4));'));
            % CUTHILL-MCKEE REORDERING
                %eval(strcat('A',num2str(k),' = A',num2str(k),'(symrcm(A',num2str(k),'),symrcm(A',num2str(k),'));'));
            % GAUSS-JORDAN ELIMINATION
                %eval(strcat('rref(A',num2str(k),');'));
            % MATLAB TRIANGLE SOLVER
                eval(strcat('L''\(L\f(1:size(A',num2str(k),')));'));
            % PCG WITH CHOLESKY FACTORIZATION
                %eval(strcat('pcg(A', num2str(k), ',f(1:size(A', num2str(k),')),10e-5,1000000,L,L'');'));
            % GMRES WITH CHOLESKY FACTORIZATION
                %eval(strcat('gmres(A', num2str(k), ',f(1:size(A', num2str(k),')),10,10e-5,1000000,L,L'');'));
            % MINRES WITH CHOLESKY FACTORIZATION
                %eval(strcat('minres(A', num2str(k), ',f(1:size(A', num2str(k),')),10e-5,1000000,L,L'');'));
        num = toc;
        % Save block time
        blocktimeelapsed(k) = num;
    end
    
    % Record time elapsed for full system
    fulltimes(z) = fulltimeelapsed;
    % Record maximum block time elapsed
    blocktimes(z) = max(blocktimeelapsed);
end

% Plot max block time as a percentage of full system time
figure()
plot(blocktimes./fulltimes*100)
hold on;
grid on;
plot(0:0.1:q,ones(1,length(0:0.1:q))*mean(blocktimes./fulltimes*100))
plot(0:0.1:q,ones(1,length(0:0.1:q))*100)
legend('Time Difference','Average Difference','Break-Even')
xlabel('Trial')
ylabel('% Time elapsed')
title(['Full Time vs Max Block Time for n=',num2str(nblocks),' Sub-Blocks'])

percenttimeelapsed = mean(blocktimes./fulltimes*100)