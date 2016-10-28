 # Read about the Van der Pol oscillator at
 #   http://en.wikipedia.org/wiki/Van_der_Pol_oscillator.

 pkg load odepkg;
 
 vopt = odeset('InitialStep', 1.0e-3, 'MaxStep', 1.0e-2, 'RelTol', 1e-3);
 [T, vsol] = ode45(@fvanderpol, [0, 100], [2; 0], vopt); 

 
 % Selecao de um trecho periodico
 [val1,idx1]=cummax(vsol(:,1));
 [val2,idx2]=cummax(vsol(:,2));
 
 y1=unique(idx1); k1=length(y1);
 y2=unique(idx2); k2=length(y2);
 
 # Isolamento de um periodo
 vsolp1=vsol(y1(length(y1)-1):y1(length(y1)),1);
 vsolp2=vsol(y2(length(y2)-1):y2(length(y2)),2);
 
 clear('vsol');
 
 vsolp1=vsolp1'; vsolp2=vsolp2';
 
# dizimacao, garantir pelo menos uma amostra por perido
 diz=20;
 if diz >= min(length(vsolp1),length(vsolp2))
   printf("Dizimacao e maior que o comprimento do vetor");
   break;
 endif;
 
# define o tamanho da sequencia desejada
lensequ=10^7;

# intervalos de quantizacao
quant=[-1 1]; % niveis de quantizacao
 
# espansao das sequencias - tamanho de 10^6 (deve ser suficiente pela periodicidade das sequencias)
# os vetores q1vsol e q2vsol correspondem as sequencias para analise
steps=ceil(log2(lensequ*diz/length(vsolp1)));
for i=1:steps
  vsolp1=[vsolp1 vsolp1];
end

# quantizacao1, sequencia de saida q1vsol
vsol1=vsolp1(1:diz:length(vsolp1));
clear('vsolp1');
k=length(vsol1);
q1vsol=zeros(k,1);
idx0=vsol1<quant(1);
idx1=vsol1<quant(2);
clear('vsol1');
idx1=xor(idx0,idx1);
idx2=!(idx0 | idx1);
q1vsol(idx1)=1;
q1vsol(idx2)=2;

# espansao das sequencias - tamanho de 10^6 (deve ser suficiente pela periodicidade das sequencias)
# o vetor q1vsol corresponde a uma das sequencias para analise
steps=ceil(log2(lensequ*diz/length(vsolp2)));
for i=1:steps
  vsolp2=[vsolp2 vsolp2];
end

# quantizacao2, sequencia de saida q2vsol
vsol2=vsolp2(1:diz:length(vsolp2));
clear('vsolp2');
k=length(vsol2);
q2vsol=zeros(k,1);
idx0=vsol2<quant(1);
idx1=vsol2<quant(2);
clear('vsolp2');
idx1=xor(idx0,idx1);
idx2=!(idx0 | idx1);
q2vsol(idx1)=1;
q2vsol(idx2)=2;
