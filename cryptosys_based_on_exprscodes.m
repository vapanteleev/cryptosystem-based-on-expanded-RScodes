clc;
clearvars;
p=5; %field characteristic
m=3; %extension power
%k/n>(1-λ/m)
k=30; %message kength
n=p^m-1; %codeword length
t=(n-k)/2;
prim_poly=[2 3 0 1]; %primitive polynomial




% vector_representation = gftuple((0:p^m-2)',prim_poly,p);

                

% ________________________________ CREATING G ___________________________________________________________________

% G={}; %cell array
% for i =1:k-1
%     for j =1:n-1
%         G{i,j}=mod(accompanying(prim_poly)^i*j,p);
%     end
%     
%     
% end
% G=cell2mat(G);
%  A = ones((k-1)*m,1);
%  B=ones(1,(n-1)*m);
%  G=horzcat(A,G);
% 
%     G=vertcat(B,G);

% G={};
% for i =1:k-1
%     for j =1:n-1
%         G{i,j}=mod(accompanying(prim_poly)^i*j,p);
%     end
%     
%     
% end
%  G=cell2mat(G);
%  A = ones((k-1)*m,1);
%  B=ones(1,(n-1)*m+1);
%  G=horzcat(A,G);
%  G=vertcat(B,G);


% _________________________________CREATING H__________________________________________________________________

H={}; %cell array of powers of accompanying matrixes of prim poly
for i = 1:(n-k) %rows of H            
    for j = 1:n-1 %columns of H
        H{i,j}=mod(accompanying(prim_poly)^i*j,p); 
    end
    
    
end
H=cell2mat(H); 
D=ones(m*(n-k),1);
 H=horzcat(D,H);

C=ones(1,m*2*(n-k)+1);

%_________________________________CHECK WETHER H ORTHOGONAL G__
% GH=mod(G*H',5);


%_________________________________________________________________________________________________________________

% H_except_once = cell2mat(struct2cell(load('Ha.mat')));
% A = ones(339,1);
% H=[A H_except_once];
% B=ones(k,1);
% C=ones(1,n-1);
% _________________________________ SHORTENING H __________________________________________________________________

%     H_short = cell2mat(struct2cell(load('Hs.mat')));

    rp=randperm(size(H,2));
    Hs=H;
    Hs(:,rp(1:n-2))=[];
    H_shorted=Hs;
%    H_short=shorteningH(H);
% ___________________________________________________________________________________________________



% T = create_T;
% T=cell2mat(T);


% ___________________________________________________________________________________________________

% ___________________________________________________________________________________________________

% _________________________KEY GENERATION____________________________________

        T={};
        T_itie_array={};
        for i=1:n
            T_itie_array{i}=randi([0, 5], 2);
        end
        Permut_T_itie_array = T_itie_array(randperm(numel(T_itie_array))) ;


    for i=1:n
            T{i,i} = T_itie_array{i};

    end

    
    for f=1:n 
        for j=1:n
            if isempty(T{f,j})
                T{f,j}=[0,0;0,0];
            end
        end
    end

    



P_sigma={};

for i=1:n
    
            P_sigma{i,i} = Permut_T_itie_array{i};

end

    
    for kk=1:n 
        for j=1:n
            if isempty(P_sigma{kk,j})
                P_sigma{kk,j}=[0,0;0,0];
            end
        end
    end
T=cell2mat(T);
P_sigma=cell2mat(P_sigma);
Q = mod(T*P_sigma,5);
H_hatch=mod(H_shorted*Q,5);

% ___________________________________________________________________________________________________

% _______________________________________ENCRYPTION__________________________________________________


lambda = 2;
prim_poly_lambda=gfprimdf(lambda,p);
% vector_representation_lambda = gftuple((0:p^lambda-2)',prim_poly_lambda,p);
sym1 = accompanying(prim_poly_lambda);
% ненулевые значения сообщения "у" в n-k позициях
y = {};
for f=1:t
     y{f}=accompanying(prim_poly_lambda);

end
for g = (t+1):n
    y{g}=[0,0;0,0];
end

    
    

y=cell2mat(y);


c=H_hatch*y';
c=mod(c,5);

%______________________DECRYPTION____________________________________________________________________

% ___________________________________________________________________________________________________

% 
% function H_shorted = shorteningH(H,n)
%     rp=randperm(size(H,2));
%     Hs=H;
%     Hs(:,rp(1:n))=[];
%    H_shorted=Hs;
%     return
%     
% end





% ___________________________________________________________________________________________________
function D = accompanying(h)
D = diag(ones(1,length(h)),1);
D(end,:) = [1 -h];
D(:,1)=[];
D(:,end)=[];
D(end-1,:)=[];
D(1,:)=[];
% D(end,:)=mod(D(end,:),13);
 D=mod(D,5);
% D = num2cell(D);

 return

 end
